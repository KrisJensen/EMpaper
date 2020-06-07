#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 11:33:09 2018
@author: loaner
"""
import pandas as pd
import networkx as nx
from pymaid import *
import copy
import numpy as np
import sys
from buildCX import buildCX
sys.setrecursionlimit(1500)

class treeNeuron():
    '''need to take a neuron, start at the root, walk upwards (keeping track of all branches)
    while adding distances
    need to have a networkx structure
    '''
    
    
    def __init__(self, neuron, inputs, outputs):
        '''give list of inputs, list of outputs. Currently only to other 'known' neurons
        can easily extend to all synapses and simply get list from catmaid'''
        
        #self.g0 = neuron.graph
        self.neuron = neuron
        self.g0 = neuron2nx(neuron)
        self.inputs = inputs
        self.outputs = outputs
        self.graph_root = neuron.nodes[ neuron.nodes.type == 'root' ].treenode_id.values[0]
        if not neuron.rootnode == 'None':
            self.rootnode = neuron.rootnode
        else:
            self.rootnode = neuron.soma
            
        if type(self.rootnode)==list:
            newroots = []
            for root in self.rootnode:
                if type(root) == tuple:
                    for subroot in root: newroots.append(subroot)
                else:
                    newroots.append(root)
            self.rootnode = newroots  #make list of rootnodes even if a subcompartment has multiple
        #rootnode has predecessors but no successors. Go to predecessors
        #self.diameter = lambda x: max(20, (500-(500-20)*x/100000))/1000 #gives diameter in um vs. dist from rootnode with mincap at 20 nm 
        self.graph = nx.DiGraph()
        self.graph.add_node(self.graph_root, Type = 'graph_root', neuron = self.neuron.neuron_name, dist = '0.0',\
                            diameter = str(self.diameter( self.neuron, self.graph_root, self.rootnode)), loc = np.array([0, 0, 0]) )
        self.rootnode = self.rootnode
        self.construct_tree(self.graph_root)
        self.synapses = []
        #walk through tree
        #if branch or synapse: add to graph
        
    def diameter(self, neuron, node, rootnode):
        '''finds our diameter at this point by considering whether or not we are axonic and then applying
        a fitted exponential decay function'''
        
        d_axon = 0.4420 #um
        d_asymptotic = 0.0250 #um
        decay = 0.000040589 #nm^-1 (note inconsistent units due to differences between pymaid & NEURON)
        if type(rootnode) == int: #we are only considering a single neuropil
            dist = dist_between(neuron, node, rootnode)
            return (d_axon-d_asymptotic)*np.exp(-decay*dist)+d_asymptotic
        else:
            newroots = []
            for root in rootnode:
                if type(root) == tuple:
                    for subroot in root: newroots.append(subroot)
                else:
                    newroots.append(root)
            rootnode
            #if we are between rootnodes, we are axonic;
            #i.e. if length from node to any two roots is less than dist  between these two
            rootdists = []
            for root in rootnode:
                rootdists.append(dist_between(neuron, node, root))
            minroot = rootnode[rootdists.index(min(rootdists))] #nearest root
            for root in rootnode:
                if not root == minroot:
                    if dist_between(neuron, node, root) < dist_between(neuron, minroot, root):
                        #if distance to another root is smaller than interoot dist, we are axonic
                        return d_axon
                    else:
                        #if interroot dist is smallest, we are dendritic
                        dist = dist_between(neuron, node, minroot)
                        return (d_axon-d_asymptotic)*np.exp(-decay*dist)+d_asymptotic #fitted decay

    def addbytype(self, dic, new, newlist, Type):
        for key in dic.keys():
            if key in new:
                if Type == 'pre':
                    dic[key] = dic[key] + [con[0] for con in newlist]
                else:
                    dic[key] = dic[key] + [con[1] for con in newlist]
        return dic
        
    def add_synapses(self, vol):
        '''
        get locations of synapses
        '''
        cons_in = {\
        'PEN1' : [],\
        'PEN2' : [],\
        'EPG' : [],\
        'PEG' : [],\
        'D7' : [],\
        'R' : []}
        cons_out = copy.deepcopy(cons_in)
        in_neuron = {}
        out_neuron = {}
        for pre, pre_dic in vol.cn_tables.items():
            if pre == self.neuron.neuron_name:
                print('adding pre synapses')
                for post, syns in pre_dic.items():
                    cons_out = self.addbytype(cons_out, post, syns, 'pre')
                    out_neuron[post] = [con[0] for con in syns]
                    if post == self.neuron.neuron_name:
                        print('adding self synapses')
                        cons_in = self.addbytype(cons_in, pre, pre_dic[self.neuron.neuron_name], 'post')
                        in_neuron[pre] = [con[1] for con in pre_dic[self.neuron.neuron_name]]
            else:
                print('adding post synapses')
                cons_in = self.addbytype(cons_in, pre, pre_dic[self.neuron.neuron_name], 'post')
                in_neuron[pre] = [con[1] for con in pre_dic[self.neuron.neuron_name]]
                
        self.out_bytype = cons_out
        self.in_bytype = cons_in
        self.out_byneuron = out_neuron 
        self.in_byneuron = in_neuron
        
    def add_node(self, nc, nm1, w, Type):
        '''add node to graph'''
        d = self.diameter( self.neuron, nc, self.rootnode )
        vec = np.random.rand(3) #move in one qudrant only to build a tree
        loc = self.graph.nodes[nm1]['loc']+vec/np.linalg.norm(vec)*w/1000 #new location is old location plus random vector of length w in um (distance from nm1)
        if not nc in self.graph.nodes:
            #dist is distance from parent
            self.graph.add_node(nc, Type = Type, neuron = self.neuron.neuron_name, dist=str(np.round(w/1000,5)), diameter = str(np.round(d,4)), loc = loc)
        else: print('node ', nc, ' already in graph')
        if not nc == nm1:   
            self.graph.add_edge(nm1, nc, dist=str(np.round(w,2)), diameter = str(np.round(d,4)), Type = 'electrical') #make neuron bidirectional
            self.graph.add_edge(nc, nm1, dist=str(np.round(w,2)), diameter = str(np.round(d,4)), Type = 'electrical')
        
        return nc, 0.0
        
    def construct_tree(self, tid, child = 'None', inputs = 'default', outputs = 'default', init = True):
        '''use the child argument to specify which child to follow if starting at branchpoint'''
        w = 0.0
        nc = tid
        nm1 = tid #last node added to self.graph  
        if init:
            outputs = copy.deepcopy(self.outputs)
            inputs = copy.deepcopy(self.inputs)
        if child == 'None': ps = list( self.g0.predecessors(tid) )
        else: ps = [child] 
        while len(ps) == 1: #only 1 child, add this to graph  
            w += self.g0.get_edge_data(ps[0], nc)['weight'] #update geodesic distance
            nc = ps[0] #update current neuron
            ps = list( self.g0.predecessors(nc) )#get new children
            if nc in inputs:#figure out how to check this quickly
                nm1, w = self.add_node(nc, nm1, w, 'input') #add node, update nm1 to nc and set w=0. add in/out later    
                inputs.pop(inputs.index(nc))
                print(self.neuron.neuron_name+' added input: '+str( ( 1 - len(inputs) / len(self.inputs) ) * 100 )+'%')
            elif nc in outputs:#figure out how to check this quickly
                nm1, w = self.add_node(nc, nm1, w, 'output')
                outputs.pop(outputs.index(nc))
            elif nc == self.rootnode:
                nm1, w = self.add_node(nc, nm1, w, 'rootnode')
            else:
                try:
                    if nc in self.rootnode: #if multiple rootnodes
                        nm1, w = self.add_node(nc, nm1, w, 'rootnode')
                except TypeError: #single rootnode
                    if nc == self.rootnode: #if multiple rootnodes
                        nm1, w = self.add_node(nc, nm1, w, 'rootnode')
                
        if len(ps) == 0:#we've reached a terminal node and can stop building  
            nm1, w = self.add_node(nc, nm1, w, 'terminal') 
            print('end of branch at ', nc)
            return
        if len(ps) > 1:
            nm1, w = self.add_node(nc, nm1, w, 'branch') #add node, update nm1 to nc and set w=0 
            for child in ps:
                self.construct_tree(nc, child = child, inputs=inputs, outputs=outputs, init=False)
            print('finished branchpoint at ', nc)
            return
    
    def write_line_neuron(self, n, node, parent, branch, f):
        '''
        write a single line of the output file
        n and parent are line numbers, node and branch are tids
        '''
        diameter = (float(self.graph.nodes[node]['diameter'])+float(self.graph.nodes[branch]['diameter']))/2 #segment width is avg
        #dist is distance from parent
        f.write('\n'+str(n)+' '+self.graph.nodes[node]['dist']+' '+str(diameter)+' '+str(parent)+' '+str(node))
        
        if node == self.rootnode: self.neuron_nroot = n
        return n+1
    
    def write_neuron(self, filename = 'tree'):
        with open(filename+'.neuron', 'w') as f:
            f.write('%neuron '+self.neuron.neuron_name) #write some preliminary info
            f.write('\n#User ktj21')
            f.write('\n#Jayaraman lab')
            #number, length, diameter, parent
            f.write('\n0 0.001 '+self.graph.nodes[self.graph_root]['diameter']+' -1 '+str(self.graph_root)) #make rootnode a 'point' section
            added = [self.graph_root]
            n = self.new_branch_neuron(self.graph_root, 1, 0, f, added)
            
            if type(self.rootnode) == int:
                f.write('\n%root ' + str(self.rootnode))
            else:
                f.write('\n%root')
                for root in self.rootnode:
                    f.write(' '+str(root))
                
            for Type, cons in self.out_bytype.items():
                f.write('\n%pre ' + Type)
                for con in cons: f.write(' '+str(con))
                
            for Type, cons in self.in_bytype.items():
                f.write('\n%post ' + Type)
                for con in cons: f.write(' '+str(con))
                
            for Type, cons in self.out_byneuron.items():
                f.write('\n%pre ' + Type)
                for con in cons: f.write(' '+str(con))
                
            for Type, cons in self.in_byneuron.items():
                f.write('\n%post ' + Type)
                for con in cons: f.write(' '+str(con))
        return
    
    def new_branch_neuron(self, branch, n, nparent, f, added):
        '''construct new branch. n counts segment number. branch'''
        for suc in self.graph.successors(branch):
            if not suc in added:
                added.append(suc)
                n = self.write_line_neuron(n, suc, nparent, branch, f) #n is line number, suc is actual node id, nparent is line number of parent, branch is node id of parent
                print(suc, n, n-1)
                n = self.new_branch_neuron(suc, n, n-1, f, added)
        return n

    
    def get_str(self, loc):
        return str(loc[0])+' '+str(loc[1])+' '+str(loc[2])+' '
    
if __name__ == '__main__':
    
    mydir = '../NEURON/neurons/' #directory to store these in
    
    show = True #show the neurons as we import them
    vol = buildCX(show=show, from_pickled = True) #construct a 'volume' instance from the neurons
    
    for neuron in ['PEN2-5R', 'PEN2-6Ra', 'PEN2-6Rb', 'PEN2-7R',
                  'PEN2-5L', 'PEN2-6La', 'PEN2-6Lb', 'PEN2-7L']:
        
        vol.get_cons_by_neuron()   
        n = copy.copy( vol.get_neuron(neuron) )
        t = treeNeuron(n, vol.neuron_inputs[n.neuron_name], vol.neuron_outputs[n.neuron_name]) #construct a tree
        t.add_synapses(vol) #get synapses from neuropil
        t.write_neuron(mydir+n.neuron_name+'_all') #write to file
    
