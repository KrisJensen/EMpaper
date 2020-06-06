#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:36:25 2019

@author: kris


code for making Dan's diagram
"""

from modelCX import modelCX
import pickle
import numpy as np
from buildCX import buildCX, get_connection_distributions
from estimate_synapses_distribution import get_connection_estimate_EB_wedge, get_connection_estimate_PB


def get_D7_D7_cons(vol):
    PB = vol.subvolumes['PB']
    
    groups = [['D7-5R4La', 'D7-5R4Lb'], ['D7-6R3La', 'D7-6R3Lb'], ['D7-7R2L'], ['D7-2R7L'], ['D7-1R9R8L']]
    
    concounts = []
    
    for i in range(len(groups)):
        g1 = groups[i]
        for j in range(i+1, len(groups)):
            g2 = groups[j]
            for n1 in g1:
                for n2 in g2:
                    try: concounts.append(len(PB.cn_tables[n1][n2]))
                    except KeyError: concounts.append(0)
                    try: concounts.append(len(PB.cn_tables[n2][n1]))
                    except KeyError: concounts.append(0)
                    
                    print(n1, n2, concounts[-2], concounts[-1])
    
    print('mean:', np.mean(concounts))
    return concounts

def write_combined_cons(continuous = False, fname = 'extrapolated/PB_EB_extrapolated', shift = 0, regular_D7s = False):
    groups = ['EPG', 'PEN1', 'D7']
    
    mod = modelCX(groups = groups, shift = shift, regular_D7s = regular_D7s)
    fname += str(shift)
    if regular_D7s: fname+= '_reg'
    
    consEB = pickle.load( open('extrapolated/avg_connections_EB_pickled', 'rb') )
    
    wsPEN = { '2R':5, '2L':5, '3R':6, '3L':4, '4R':7, '4L':3, '5R':8, '5L':2,\
             '6R':1, '6L':1, '7R':2, '7L':8, '8R':3, '8L':7, '9R':4, '9L':6 }
    
    wsEPG = { '1R':5, '1L':5, '2R':6, '2L':4, '3R':7, '3L':3, '4R':8, '4L':2,\
             '5R':1, '5L':1, '6R':2, '6L':8, '7R':3, '7L':7, '8R':4, '8L':6 }
    
    wsPENw = { '2R':(51,52), '2L':(51,52), '3R':(61,62), '3L':(41,42), '4R':(71,72), '4L':(31,32), '5R':(81,82), '5L':(21,22),\
             '6R':(11,12), '6L':(11,12), '7R':(21,22), '7L':(81,82), '8R':(31,32), '8L':(71,72), '9R':(41,42), '9L':(61,62) }
    
    wsPEGw = { '1R':(51,52), '1L':(51,52), '2R':(61,62), '2L':(41,42), '3R':(71,72), '3L':(31,32), '4R':(81,82), '4L':(21,22),\
             '5R':(11,12), '5L':(11,12), '6R':(21,22), '6L':(81,82), '7R':(31,32), '7L':(71,72), '8R':(41,42), '8L':(61,62) }
    
    wsEPGw = { '1R':(52,), '1L':(51,), '2R':(62,), '2L':(41,), '3R':(72,), '3L':(31,), '4R':(82,), '4L':(21,),\
             '5R':(12,), '5L':(11,), '6R':(22,), '6L':(81,), '7R':(32,), '7L':(71,), '8R':(42,), '8L':(61,) }

    for g1 in groups: #don't need D7 for EB
        for g2 in groups:
            for n1 in mod.groups[g1]:
                for n2 in mod.groups[g2]:
                    if n1!=n2:
                        try:
                            
                            if 'R' in g1 and 'R' in g2:
                                continue
                            
                            elif 'R' in g2:
                                mod.add_connection(n1, n2, consEB[g1]['output'][g2])
                            
                            elif 'R' in g1:
                                mod.add_connection(n1, n2, consEB[g2]['input'][g1])
                                
                            else:
                                ws = []
                                for n in (n1, n2): #find wedges for n1 and n2
                                    if 'PEN' in  n: ws.append( wsPENw[n[5:7]] )
                                    elif 'GE' in n:
                                        if n[3] == 'r': ws.append( ( int(n[4]+'1'), ) ) #if it goes to right Gall it overlaps with EPG-Lx (left pB goes to right GAll)
                                        elif n[3] == 'l': ws.append( ( int(n[4]+'2'), ) )
                                        else: print('GE not recognized', n)
                                    elif 'PEG' in n: ws.append( wsPEGw[n[4:6]] )
                                    elif 'EPG' in n: ws.append( wsEPGw[n[4:6]] )
                                    else:
                                        {}[0] #force key error and zero connectivity
                                    
                                same = False
                                neighbor = False
                                
                                for w1 in ws[0]: #check if there is overlap
                                    for w2 in ws[1]: 
                                        if w1 == w2: same = True #there is overlap
                                        
                                        elif str(w1)[0] == str(w2)[0]: neighbor = True #same tile; neighbors
                                        
                                        elif (int(str(w2)[0])-int(str(w1)[0]) == 1 or int(str(w1)[0])-int(str(w2)[0]) == 7)\
                                            and int(str(w1)[1]) - int(str(w2)[1]) == 1: #e.g. 52 and 61, 82 and 11 are neighbours
                                            neighbor = True
                                            
                                        elif (int(str(w1)[0])-int(str(w2)[0]) == 1 or int(str(w2)[0])-int(str(w1)[0]) == 7)\
                                            and int(str(w2)[1]) - int(str(w1)[1]) == 1: #e.g. 61 and 52, 11 and 82 are neighbours
                                            neighbor = True
                                        
                                if same:
                                    mod.add_connection(n1, n2, consEB[g1]['output'][g2+'_s']) #if there us any overlap, we consider the fragments to the in the same wedge
                                    
                                elif neighbor: #neighbor wedges
                                    mod.add_connection(n1, n2, consEB[g1]['output'][g2+'_n']) #if there is any NN interactions, we consider them to be neighboring
                                    
                        except KeyError:
                            #print('keyerror', n1, n2)
                            mod.add_connection(n1, n2, 0)
                                
    

    consPB = pickle.load( open('extrapolated/avg_connections_PB_pickled', 'rb') )
        
    for g1 in groups:
        for g2 in groups:
            for n1 in mod.groups[g1]:
                for n2 in mod.groups[g2]:
                    if n1!=n2:
                        try:

                            gs = []
                            for n in (n1, n2): #find glomeruli for n1 and n2
                                if 'PEN' in  n: gs.append( n[5:7] )
                                
                                elif 'D7' in n:
                                    if '9R' in n or '9L' in n: gs.append( (n[-3:-1],n[-5:-3],n[-7:-5]) )
                                    else: gs.append( (n[-3:-1],n[-5:-3]) )
                                
                                else: gs.append( n[4:6] )
                                
                            if type(gs[0]) == tuple and type(gs[1]) == tuple: #D7 onto D7
                                mod.connections[n1][n2] += 5.5 #just take average of D7-D7 connections from function above
        
                            if gs[0] == gs[1]: #same glomerulus; should we add a connection from EPG-1L t PEN1-9L?
                                mod.connections[n1][n2] += consPB[g1]['output'][g2]
                                
                            elif (continuous and ( (gs[0] == '1L' and gs[1] == '9L') or (gs[0] == '1R' and gs[1] == '9R') ) ): #add EPG to PEN connection to complete loop
                                mod.connections[n1][n2] += consPB[g1]['output'][g2]
                                    
                            elif type(gs[0]) == tuple: #n1 is d7, outputs to 2 glomeruli
                                if gs[1] in gs[0]: #n2 gets input from D7
                                    mod.connections[n1][n2] += consPB[g2]['input'][g1]
                            elif type(gs[1]) == tuple: #n2 id D7, gets input from remaining glomeruli
                                if not gs[0] in gs[1]:
                                    mod.connections[n1][n2] += consPB[g1]['output'][g2]
                                           
                        except KeyError:
                            mod.connections[n1][n2] += 0
                            
                                
    mod.connectivity(groups1='all', groups2='all', file='connectivity_matrices/'+fname,\
                       heatmap=fname, grouping='default', size = (25, 25) )
    
    
    return mod, consEB, consPB
    
if __name__ == '__main__':

    show = True #show the neurons as we import them
    vol = buildCX(show=show, from_pickled = True) #construct a 'volume' instance from the neurons

    vol.subvolumes['EB'].get_connection_distribution(relation = 'input')
    vol.subvolumes['EB'].get_connection_distribution(relation = 'output')
    est, avgCon = get_connection_estimate_EB_wedge(vol.subvolumes['EB'])
    pickle.dump(avgCon, open('extrapolated/avg_connections_EB_pickled', 'wb') )
    
    vol.subvolumes['PB'].get_connection_distribution(relation = 'input')
    vol.subvolumes['PB'].get_connection_distribution(relation = 'output') 
    est, avgCon = get_connection_estimate_PB(vol.subvolumes['PB'])
    pickle.dump(avgCon, open('extrapolated/avg_connections_PB_pickled', 'wb') )
    
    mod, consEB, consPB = write_combined_cons(continuous = True, fname = 'extrapolated/PB_EB_cont_extrapolated')

    
    

