#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 13:32:56 2018

@author: loaner

"""


import pymaid
import matplotlib.pyplot as plt
import numpy as np


def get_dist(nodes, tid1, tid2):
    try: return int(np.sqrt( (nodes.loc[tid1].x-nodes.loc[tid2].x)**2 +\
                              (nodes.loc[tid1].y-nodes.loc[tid2].y)**2 +\
                              (nodes.loc[tid1].z-nodes.loc[tid2].z)**2 ))
    except TypeError:
        print('node does not exist')
        return 0
    
def remove_node(nodes, tid):
    cids = nodes[ nodes.parent_id == tid ].index.values
    try:
        pid = nodes.loc[tid, 'parent_id']
        nodes = nodes.drop(tid)
        for cid in cids:
            nodes.loc[cid, 'parent_id'] = pid        
        return nodes
    except TypeError:
        print('removing terminal node')
        print(nodes.loc[tid])
        nodes = nodes.drop(tid)
        for cid in cids:
            nodes.loc[cid, 'parent_id'] = None
        return nodes

def dist_despike(nodes, threshold, show):
    
    repeat = False
    dists = []
    for tid in nodes.index.values:
        try: pid = nodes.loc[tid].parent_id    
        except KeyError: continue    
        if not pid == None:        
            dist = int(np.sqrt( (nodes.loc[tid].x-nodes.loc[pid].x)**2 +\
                              (nodes.loc[tid].y-nodes.loc[pid].y)**2 +\
                              (nodes.loc[tid].z-nodes.loc[pid].z)**2 ))        
            if dist > threshold:            
                print('\ntid', tid, 'pid', pid, 'dist', dist)
                n = 0
                m = 0
                remove = []        
                tid0 = tid
                pid0 = pid            
                print()
                try:
                    cid = nodes[ nodes.parent_id == tid ].index.values[0]            
                    while get_dist(nodes, tid, cid) < threshold and n < 4:
                        n += 1
                        remove.append(tid)
                        tid = cid
                        try: cid = nodes[ nodes.parent_id == tid ].index.values[0]
                        except:
                            print('reached end, removing nodes')
                            n = 4
                    if len(remove) == 4:
                        print('not removing children')
                        remove = []
                    elif len(remove) == 0:
                        print('removing node', tid0)
                        nodes = remove_node(nodes, tid0)
                        repeat = True
                    else:
                        print('removing children', remove)
                        repeat = True
                        for tid in remove:
                            nodes = remove_node(nodes, tid)
                except IndexError:
                    print('no children, remove node')
                    nodes = remove_node(nodes, tid)
                    repeat = True
                try:
                    remove = []
                    ppid = nodes.loc[pid0].parent_id
                    while get_dist(nodes, pid, ppid) < threshold and m < 6:
                        m += 1
                        remove.append(pid)
                        pid = ppid
                        try: ppid = nodes.loc[pid].parent_id
                        except:
                            print('reached end, removing nodes')
                            m = 6
                    if m == 6:
                        print('not removing parents')
                        remove = []
                    elif m == 0:
                        print('removing parent', pid0)
                        nodes = remove_node(nodes, pid0)
                        repeat = True
                    else:
                        print('removing parents', remove)
                        for pid in remove:
                            nodes = remove_node(nodes, pid)
                        repeat = True
                except TypeError:
                    print('no grandparent, remove parent')
                    nodes = remove_node(nodes, pid0)
                    repeat = True
            else: dists.append(dist)
    if show:        
        plt.hist(dists)
        plt.show()
    if not repeat:
        print('no change, moving on')
        repeat = False
    return nodes, repeat

def despike(neuron, sigma = 4, threshold = 2000, show = True):
    if not sigma == 'None': pymaid.despike_neuron(neuron, sigma, inplace = True, max_spike_length = 5 )
    repeat = True
    if show:        
        neuron.plot2d()
        plt.show()
    nodes = neuron.nodes.set_index('treenode_id')
    while repeat:
        print('new round of despiking')
        nodes, repeat = dist_despike(nodes, threshold, show)
        neuron.nodes = nodes
        # Reassign treenode table
        neuron.nodes = nodes.reset_index(drop=False)
        # The weights in the graph have changed, we need to update that
        neuron._clear_temp_attr()
        if show:
            neuron.plot2d()
            plt.show()
    return neuron

    
    