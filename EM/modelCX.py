#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 16:15:03 2018

@author: loaner

Write script for 'completing' the central complex by extrapolating our average neurons.

"""

from volume import volume
import matplotlib.pyplot as plt
import numpy as np
import pickle
import networkx as nx


class modelCX():
    
    
    def __init__(self, groups='all', shift = 0, regular_D7s = False):
        ''' shift determines which PEN neuron to start wth (default PEN1-2R)
        if regular_D7s, leave out the Delta7s that don't fit the wiring pattern'''
    
        neurons = {}
        
        neurons['PEN1'] = ['PEN1-2R', 'PEN1-3R', 'PEN1-4R', 'PEN1-5R', 'PEN1-6R', 'PEN1-7R', 'PEN1-8R', 'PEN1-9R',\
            'PEN1-9L', 'PEN1-8L', 'PEN1-7L', 'PEN1-6L', 'PEN1-5L', 'PEN1-4L', 'PEN1-3L', 'PEN1-2L']
        
        for s in range(shift): neurons['PEN1'] = [neurons['PEN1'].pop()]+neurons['PEN1']
    
        neurons['PEN2'] = [ PEN[0:3]+'2'+PEN[4:] for PEN in neurons['PEN1']]
    
    
        neurons['PEG'] = ['PEG-1R', 'PEG-1L', 'PEG-2R', 'PEG-2L', 'PEG-3R', 'PEG-3L', 'PEG-4R', 'PEG-4L', 'PEG-5R', 'PEG-5L',\
             'PEG-6R', 'PEG-6L', 'PEG-7R', 'PEG-7L', 'PEG-8R', 'PEG-8L']
        
    
        neurons['EPG'] = ['EPG-1La','EPG-1Lb','EPG-1Lc', 'EPG-1Ra','EPG-1Rb','EPG-1Rc', 'EPG-8La','EPG-8Lb','EPG-8Lc',\
            'EPG-2Ra','EPG-2Rb','EPG-2Rc', 'EPG-7La','EPG-7Lb','EPG-7Lc', 'EPG-3Ra','EPG-3Rb','EPG-3Rc',\
            'EPG-6La','EPG-6Lb','EPG-6Lc', 'EPG-4Ra','EPG-4Rb','EPG-4Rc', 'EPG-5La','EPG-5Lb','EPG-5Lc',\
            'EPG-5Ra','EPG-5Rb','EPG-5Rc', 'EPG-4La','EPG-4Lb','EPG-4Lc', 'EPG-6Ra','EPG-6Rb','EPG-6Rc',\
            'EPG-3La','EPG-3Lb','EPG-3Lc', 'EPG-7Ra','EPG-7Rb','EPG-7Rc', 'EPG-2La','EPG-2Lb','EPG-2Lc',\
            'EPG-8Ra','EPG-8Rb','EPG-8Rc']
        
        neurons['GE'] = ['GE-r1', 'GE-l1', 'GE-r2', 'GE-l2','GE-r3', 'GE-l3', 'GE-r4', 'GE-l4',\
           'GE-r5', 'GE-l5', 'GE-r6', 'GE-l6','GE-r7', 'GE-l7', 'GE-r8', 'GE-l8']
    
        neurons['NO-LAL-G'] = ['NO-LAL-G-ra', 'NO-LAL-G-rb', 'NO-LAL-G-la', 'NO-LAL-G-lb']
    
        neurons['R'] = ['R']
        
        neurons['D7'] = ['D7-9R1R8L', 'D7-8R1L', 'D7-7R2L', 'D7-6R3L', 'D7-5R4L',\
                       'D7-4R5L', 'D7-3R6L', 'D7-2R7L', 'D7-1R8L', 'D7-8R1L9L']
        
        if regular_D7s:
            neurons['D7'] = ['D7-7R2L', 'D7-6R3L', 'D7-5R4L', 'D7-4R5L', 'D7-3R6L', 'D7-2R7L']
        
        neurons['D7'].reverse()
        labels = ['a', 'b', 'c', 'd']
        neurons['D7'] = list(np.concatenate([ [n+lab for lab in labels] for n in neurons['D7'] ]))
        
        self.neurons = neurons
        
        self.groups = {}
        
        if groups == 'all':
            groups = ['EPG', 'PEG', 'PEN1', 'PEN2', 'NO-LAL-G', 'D7', 'R', 'GE']
        
        for g in groups:
            self.groups[g] = neurons[g]
            
        self.subvolumes = {}
        
        self.connections = {}
        
        print(groups)
        for g in groups:
            for n in neurons[g]:
                self.connections[n] = {}
                for g2 in groups:
                    #if g2 == 'PEN2': print('adding PEN2')
                    for n2 in neurons[g2]:

                        self.connections[n][n2] = 0
                        
        #print(self.connections['PEN1-2R'].keys())
        
    def add_subvolume(self, name, groups):
        self.subvolumes[name] = modelCX(groups)
        
    def add_connection(self, pre, post, count):
        if str(count) == 'nan': count = 0
        #if not ('GE' in pre or 'NO' in pre): print('added', pre, post, count)
        self.connections[pre][post] = count
        
    def connectivity(self, groups1='all', groups2='all', file='None',\
                       heatmap='None', grouping='default', size = (15, 15), neurons = 'all'):
        '''given two groups of neurons, give matrix of all connections FROM group1 TO group2
        file specifies name of file to write data to, method whether we want to create new file or append it
        heatmap gives name of file to write heatmap to. From connectors specifies whether or not we're using a frgment
        or an intact neuron'''
        
        if groups1 == 'all': groups1 = self.groups.keys()
        if groups2 == 'all': groups2 = self.groups.keys()
        group1 = []
        group2 = []
        for g in groups1:
            group1 += self.groups[g]
        for g in groups2:
            group2 += self.groups[g]
            
        if not neurons == 'all':
            group1 = neurons
            group2 = neurons
        
        data=[]
        tot=0.0

        
        for pre in group1: #consider each presynaptic neuron                     
            
            data.append([])
            for post in group2:
                #print(pre, post, self.connections[pre][post])
                data[-1].append( self.connections[pre][post]) #find connections to post
            
        if not file == 'None': #we write our result to a file
            with open(file+'.txt', 'w') as f:

                f.write('{:<10}'.format(' '))
                for post in group2: f.write( '{:<10}'.format( post[:10] )+' ' )
                f.write( '\n' )

                for i, pre in enumerate(group1):
                    f.write( '{:<10}'.format( pre[:10] )+' ' )
                    
                    for j, post in enumerate(group2):
                        f.write( '{:<10}'.format( str(data[i][j])[:10] )+' ' )
                    f.write('\n')
                    
                pickle.dump(data, open(file+'.pickled', 'wb'))
                    
        if not heatmap == 'None': #make a heatmap
            xlabel  = [ post for post in group2 ]
            ylabel = [ pre for pre in group1 ]
            xticks = len(xlabel)
            yticks = len(ylabel)
            
            vmax = min(100, np.array(data).max())
            
            data.reverse() #plot down to up
            ylabel.reverse() #pllot down to up
            self.heatmap(data, xticks=xticks, yticks=yticks, xlabel=xlabel, ylabel=ylabel,\
                         title=heatmap, size = size, vmax = vmax)
            
        return

    def heatmap(self, data,  xticks=None, yticks=None, xlabel=None, ylabel=None,\
                title=None, vmin=None, vmax=None, size = (15,15) ):
        '''construct heatmap'''
        
        fig, ax = plt.subplots(figsize=size)
        im = ax.imshow(data, cmap='coolwarm', vmin=vmin, vmax=vmax) #generate heatmap; could add more options
        cb = fig.colorbar(im, shrink=1.1, aspect=5 ) #add colorbar
        cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=int(2*size[0]/3))

        ax.set_xticks(np.arange(xticks))
        ax.set_yticks(np.arange(yticks))
        ax.set_xticklabels(xlabel) #add cellnames
        ax.set_yticklabels(ylabel) #add cellnames
        
        ax.xaxis.set_label_position('top') 
        ax.xaxis.tick_top()

        # Rotate the tick labels and set their alignment
        plt.setp(ax.get_xticklabels(), rotation=-45, ha="right", va='bottom', rotation_mode="anchor") #make it look nice


        plt.xlabel(title, fontsize = int(size[0]) )
        fig.tight_layout()
        
        plt.savefig('connectivity_matrices/'+title+'.png', dpi=360)
        plt.savefig('connectivity_matrices/'+title+'.pdf', dpi=360)
        
        plt.show()
        return
    

    
    