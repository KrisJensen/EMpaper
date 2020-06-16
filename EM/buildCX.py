#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 09:55:28 2018

@author: loaner
"""

from volume import volume
from pymaid import *
import pickle

def buildCX(show=True):
    '''create a volume with our neurons of interest'''
    
    vol = volume()
    groups = ['PEN1', 'PEN2', 'EPG', 'PEG', 'D7', 'R']  #neuron types to consider

    #rootnodes for the different neurons and neuropils    
    rootnodes_PEN1 = [ [3126471,3593898,2185557,2181724,2210880,2857643,3513461,12603123],\
                      [5837503,5398965,1948755,2181016,(2215133,2215382),2214922,5627995,3501313],\
                      [6870055,8103838,12900279,12900633,12900783,2418633,7252018,20883674]  ] #PB EB NO
    
    rootnodes_PEN2 = [ [3497680,4178983,2178615,2163389,2207249,2244795,3509336,5237022],\
                     [5908353,3498186,(2177781,2177647),2161670,2424373,2665948,5626946,5464801],\
                     [7019593,12905027,20799546,3299024,2206510,2203193,5564596,5584414]  ] #PB EB NO
    
    
    rootnodes_EPG = [ [ 1753050, 1748660, 3449161, 2905038, 13237900, 2262262, 6902310, 3448054, 27010310, 3554640],\
                     [ 6121805, 5936707, 3466125, 2265604, 2264211, 2864637, 16815483, 3448153, 3515658, 3380387],\
                     [ 2452093, 3309827, 2500612, 2265429, 2263551, 3302643, 20874147, 5101718, 20851318, 7518649] ] #PB EB G
    
    rootnodes_PEG = [ [2102880,4144157,3090295,3982268],\
                     [1777240,4145168,3400032,19562579],\
                     [1777173,4574589,4832222,5039997] ] #PB EB G
    
    rootnodes_D7 = [ [3008985,5289720,3426438,3241000,2849576,6295467,4061178 ],\
                    [2906483,5110398,5269799,3142774,5288727,3376503,4060762],\
                    [3024153,5109132,5997836,6930815,6632029,3376805,3438887] ]
                     
    rootnodes_R = [[21262742,21298056,21477619,22701513,23145018,23151362,23211427],\
                   [21271375,21295105,21487492,22267264,23117229,23161634,23222690]] #add this
    
    rootnodes = [rootnodes_PEN1, rootnodes_PEN2, rootnodes_EPG, rootnodes_PEG, rootnodes_D7, rootnodes_R]
    
    #load locally stored neurons
    for i, group in enumerate(groups):
        print('\n\nAdding group ', group)
        neurons = pickle.load( open('pickled_neurons/'+group+'_pickled', 'rb' ) )
        for j, neuron in enumerate(neurons):
            vol.add_neuron(neuron, subgroup = group,\
            rootnode = [ rootnodes[i][k][j] for k in range(len(rootnodes[i])) ])
    D7n = pickle.load( open('pickled_neurons/D7_pickled', 'rb' ) )
            
    subvolumes = ['PB', 'EB', 'NO'] 
    for v in subvolumes: vol.subvolume(v) #constuct subvolumes

    #split out neurons according to subvolumes    
    PEN1s = []  
    PEN2s = []
    #order PB, EB, NO

    PEN1s.append(vol.split('PEN1-5R', 3126155, 3125538, show=show))
    PEN1s.append(vol.split('PEN1-5L', 3075961, 3075572, show=show))
    PEN1s.append(vol.split('PEN1-6Ra', 2185424, 2184847, show=show))
    PEN1s.append(vol.split('PEN1-6Rb', 2181545, 2180692, show=show))
    PEN1s.append(vol.split('PEN1-6La', 2210651, 2210229, show=show))
    PEN1s.append(vol.split('PEN1-6Lb', 2857292, 2212325, show=show))
    PEN1s.append(vol.split('PEN1-7R', 3513271, 3513176, show=show))
    PEN1s.append(vol.split('PEN1-7L', 12600793, 3501452, show=show))
    
    PEN2s.append(vol.split('PEN2-5R', 3497600, 3498016, show=show))
    PEN2s.append(vol.split('PEN2-5L', 3498351, 3498551, show=show))
    PEN2s.append(vol.split('PEN2-6Ra', 2177934, 2177486, show=show))
    PEN2s.append(vol.split('PEN2-6Rb', 2163182, 2162408, show=show))
    PEN2s.append(vol.split('PEN2-6La', 2206952, 2206635, show=show))
    PEN2s.append(vol.split('PEN2-6Lb', 3389770, 12791980, show=show))
    PEN2s.append(vol.split('PEN2-7R', 3508573, 3508198, show=show))
    PEN2s.append(vol.split('PEN2-7L', 4197201, 5584341, show=show))
    
    EPGs=[]
    PEGs=[]
    #Order PB, EB, G
    EPGs.append(vol.split('EPG-5Ra', 1753455, 1771173, show=show, Type = 'fj'))
    EPGs.append(vol.split('EPG-5Rb', 1732713, 1737390, show=show, Type = 'fj'))
    EPGs.append(vol.split('EPG-5Rc', 1736306, 1741834, show=show, Type = 'fj'))
    EPGs.append(vol.split('EPG-5La', 2264515, 2265192, show=show, Type = 'fj'))
    EPGs.append(vol.split('EPG-5Lb', 2263776, 2263478, show=show, Type = 'd'))
    EPGs.append(vol.split('EPG-5Lc', 2261872, 2261504, show=show, Type = 'fj', reverse = True))
    EPGs.append(vol.split('EPG-6Ra', 3437960, 1760877, show=show, Type = 'fj'))
    EPGs.append(vol.split('EPG-6Rb', 3448088, 3448328, show=show, Type = 'new'))
    EPGs.append(vol.split('EPG-6La', 3515565, 3515922, show=show, Type = 'fj'))
    EPGs.append(vol.split('EPG-6Lb', 3378068, 3516553, show=show, Type = 'fj'))

    PEGs.append(vol.split('PEG-5R', 1777475, 1777052, show=show))
    PEGs.append(vol.split('PEG-5L', 4144656, 4573004, show=show))
    PEGs.append(vol.split('PEG-6R', 3399687, 3400550, show=show))
    PEGs.append(vol.split('PEG-6L', 3983410, 3984407, show=show))
    
    #add to subvolumes (neuropils) as appropriate
    for i, PEN1 in enumerate(PEN1s):
        print('PEN1', i)
        vol.add_to_subvolume(PEN1[0], 'PB', subgroup='PEN1', rootnode = rootnodes_PEN1[0][i], roottype = 'out')
        vol.add_to_subvolume(PEN1[1], 'EB', subgroup='PEN1', rootnode = rootnodes_PEN1[1][i], roottype = 'in')
        vol.add_to_subvolume(PEN1[2], 'NO', subgroup='PEN1', rootnode = rootnodes_PEN1[2][i], roottype = 'in')
        
    for i, PEN2 in enumerate(PEN2s):
        print('PEN2', i)
        vol.add_to_subvolume(PEN2[0], 'PB', subgroup='PEN2', rootnode = rootnodes_PEN2[0][i], roottype = 'out')
        vol.add_to_subvolume(PEN2[1], 'EB', subgroup='PEN2', rootnode = rootnodes_PEN2[1][i], roottype = 'in')
        vol.add_to_subvolume(PEN2[2], 'NO', subgroup='PEN2', rootnode = rootnodes_PEN2[2][i], roottype = 'in')
    
    for i, EPG in enumerate(EPGs):
        print('EPG', i)
        vol.add_to_subvolume(EPG[0], 'PB', subgroup='EPG', rootnode = rootnodes_EPG[0][i], roottype = 'in')
        vol.add_to_subvolume(EPG[1], 'EB', subgroup='EPG', rootnode = rootnodes_EPG[1][i], roottype = 'out')
    
    for i, PEG in enumerate(PEGs):
        print('PEG', i)
        vol.add_to_subvolume(PEG[0], 'PB', subgroup='PEG', rootnode = rootnodes_PEG[0][i], roottype = 'out')
        vol.add_to_subvolume(PEG[1], 'EB', subgroup='PEG', rootnode = rootnodes_PEG[1][i], roottype = 'in')
     
    for D7 in D7n:
        vol.add_to_subvolume(D7, 'PB', subgroup = 'D7')

    for R in ['Ra','Rb','Rc','Rd','Re','Rf','Rg']:
        vol.add_to_subvolume(R, 'EB', subgroup = 'R')
    
    #load connection tables
    vol.cn_tables = pickle.load( open('pickled_neurons/all_cn_tables_LAL_R_pickled', 'rb' ) ) 
    
    #load connection tables for subvolumes (neuropils)
    for key, subvol in vol.subvolumes.items():
        subvol.cn_tables = pickle.load( open('pickled_neurons/'+key+'_cn_tables_LAL_R_pickled', 'rb' )) 
                            
    for subvolume in vol.subvolumes:
        for group in groups:
            if not group in vol.subvolumes[subvolume].groups.keys():
                vol.subvolumes[subvolume].init_subgroup(group)
    
    return vol

def connectivity(vol, show, connectors = False, table='none'):
    '''
    construct connectivity matrices for each neuropil and the full dataset
    '''
    tots = dict()
    mats = dict()
    tot, mat = vol.allconnections(file = 'matrix_all_NLG_R', heatmap='matrix_all_NLG_R', show=show, table=table, from_connectors = connectors)
    tots['all']=tot
    mats['all']=mat
    sizes = {'NO': (10, 10), 'G': (12, 12), 'PB': (16, 16), 'EB': (20, 20) } #default figure sizes
    for region in ['PB', 'EB', 'NO']: #for each neuropil
        print('\n\nnew region ', region)
        #load connections
        tot, mat = vol.subvolumes[region].allconnections(file='matrix_'+region, table=table,\
              method='w', from_connectors=connectors, heatmap='matrix_'+region, show=show, size = sizes[region])
        tots[region] = tot
        mats[region] = mat
    
    #re-plot EB with grouping according to compartment
    tot, mat = vol.subvolumes['EB'].allconnections(file='matrix_EB_grouped', size = (20, 20), table=table,\
          method='w', from_connectors=connectors, heatmap='matrix_EB_grouped', show=show, grouping = 'EB')
        
    return tots, mats

def get_connection_distributions(vol, Norm = True, include_unknown = False):
    '''
    load and plot connection distributions as histograms in each neuropil
    '''
    for relation in ['input', 'output']:
        vol.get_connection_distribution(relation = relation )
        vol.plot_connection_distribution(relation = relation, title = 'Distribution of '+relation, estimate=False, Norm = Norm, include_unknown = include_unknown)
        for name, subvolume in vol.subvolumes.items():
            subvolume.get_connection_distribution(relation = relation )
            subvolume.plot_connection_distribution(relation = relation, title = 'Distribution of '+relation+'s '+name, estimate=False, Norm = Norm, include_unknown = include_unknown)
    

if __name__ == '__main__':

    show = True
    vol = buildCX(show=show)
    

        
        
        
        