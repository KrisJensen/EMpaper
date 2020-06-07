#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 15:14:03 2018

@author: loaner
"""
import copy
from volume import volume
import numpy as np
import pickle

def get_connection_estimate_PB(vol):
    '''go from connection distribution and try to guesstimate what the total number of connections for each is
    need to specify which region we're in for this to work. Will be quite hardcoded :(
    
    PB
    Assume we have correct number of PEN1/PEN2
    Assume we have correct number of EPG/PEG in 5/6, transfer the mean of these interactions to 7
    This implies that we have full EPG/PEG connectivity and full EPG to PEN
    Assume we have 40 D7s... Assume each covers 2 glomeruli - this gives 80 sets of output - roughly four
    per glomerulus

    NO
    Assume everything interacts with everything: scale by full number for all inputs/outputs
    
    EB
    consider same and neighboring wedges
        
    #Have now included the observation that a D7 does not get input in a region where it gives outputs.
    '''
    
    connectionEstimate = {}
        
    connections_in = copy.deepcopy(vol.connectionDistribution['input'])
    connections_out = copy.deepcopy(vol.connectionDistribution['output'])
    
        
    avgCon = {'PEN1':{}, 'PEN2':{}, 'EPG':{}, 'PEG':{}}
    for k1, g1 in avgCon.items():
        g1['input'] = {}
        g1['output'] = {}
        
    countsin = []
    countsout = []
    for EPG in vol.groups['EPG']:
        countsin.append(vol.connectionDistribution['input'][EPG.neuron_name]['tot'])
        countsout.append(vol.connectionDistribution['output'][EPG.neuron_name]['tot'])
    EPG_in = np.mean(countsin)
    EPG_out = np.mean(countsout)    

    countsin = []
    countsout = []
    for PEG in vol.groups['PEG']:
        countsin.append(vol.connectionDistribution['input'][PEG.neuron_name]['tot'])
        countsout.append(vol.connectionDistribution['output'][PEG.neuron_name]['tot'])
    PEG_in = np.mean(countsin)
    PEG_out = np.mean(countsout) 

    '''
    countsin = []
    countsout = []
    for PEN1 in vol.groups['PEN1']:
        countsin.append(vol.connectionDistribution['input'][PEN1.neuron_name]['tot'])
        countsout.append(vol.connectionDistribution['output'][PEN1.neuron_name]['tot'])
    PEN1_in = np.mean(countsin)
    PEN1_out = np.mean(countsout)

    countsin = []
    countsout = []
    for PEN2 in vol.groups['PEN2']:
        countsin.append(vol.connectionDistribution['input'][PEN2.neuron_name]['tot'])
        countsout.append(vol.connectionDistribution['output'][PEN2.neuron_name]['tot'])
    PEN2_in = np.mean(countsin)
    PEN2_out = np.mean(countsout)
    '''
    
    EPGs = [ ['EPG-6Ra', 'EPG-6Rb'], ['EPG-5Ra', 'EPG-5Rb', 'EPG-5Rc'],\
            ['EPG-5La', 'EPG-5Lb', 'EPG-5Lc'], ['EPG-6La', 'EPG-6Lb'] ]   
    
    for PEN in ['PEN1', 'PEN2']: 

        countsin = []
        countsout = []
        for P in vol.groups[PEN]:
            countsin.append(vol.connectionDistribution['input'][P.neuron_name]['tot'])
            countsout.append(vol.connectionDistribution['output'][P.neuron_name]['tot'])
        PEN_in = np.mean(countsin)
        PEN_out = np.mean(countsout)

        mean_PEN_from_EPG = np.mean( [len(vol.cn_tables['EPG-5Ra'][PEN+'-5R'])/vol.connectionDistribution['input'][PEN+'-5R']['tot'],\
                                 len(vol.cn_tables['EPG-5Rb'][PEN+'-5R'])/vol.connectionDistribution['input'][PEN+'-5R']['tot'],\
                                 len(vol.cn_tables['EPG-5Rc'][PEN+'-5R'])/vol.connectionDistribution['input'][PEN+'-5R']['tot'],\
                                 len(vol.cn_tables['EPG-5La'][PEN+'-5L'])/vol.connectionDistribution['input'][PEN+'-5L']['tot'],\
                                 len(vol.cn_tables['EPG-5Lb'][PEN+'-5L'])/vol.connectionDistribution['input'][PEN+'-5L']['tot'],\
                                 len(vol.cn_tables['EPG-5Lc'][PEN+'-5L'])/vol.connectionDistribution['input'][PEN+'-5L']['tot'],\
                                 len(vol.cn_tables['EPG-6Ra'][PEN+'-6Ra'])/vol.connectionDistribution['input'][PEN+'-6Ra']['tot'],\
                                 len(vol.cn_tables['EPG-6Rb'][PEN+'-6Ra'])/vol.connectionDistribution['input'][PEN+'-6Ra']['tot'],\
                                 len(vol.cn_tables['EPG-6Ra'][PEN+'-6Rb'])/vol.connectionDistribution['input'][PEN+'-6Rb']['tot'],\
                                 len(vol.cn_tables['EPG-6Rb'][PEN+'-6Rb'])/vol.connectionDistribution['input'][PEN+'-6Rb']['tot'],\
                                 len(vol.cn_tables['EPG-6La'][PEN+'-6La'])/vol.connectionDistribution['input'][PEN+'-6La']['tot'],\
                                 len(vol.cn_tables['EPG-6Lb'][PEN+'-6La'])/vol.connectionDistribution['input'][PEN+'-6La']['tot'],\
                                 len(vol.cn_tables['EPG-6La'][PEN+'-6Lb'])/vol.connectionDistribution['input'][PEN+'-6Lb']['tot'],\
                                 len(vol.cn_tables['EPG-6Lb'][PEN+'-6Lb'])/vol.connectionDistribution['input'][PEN+'-6Lb']['tot'] ] ) 

        for P in [PEN+'-7L', PEN+'-7R']:
            connections_in[P]['EPG'] = 3*mean_PEN_from_EPG*vol.connectionDistribution['input'][P]['tot'] #~3 EPG per glomerulus        
        avgCon[PEN]['input']['EPG'] = mean_PEN_from_EPG*PEN_in
            
        #D7 output is locaized so we take the mean value for pairs of neurons where synapses colocalize
        mean_PEN_from_D7 = np.mean( [len(vol.cn_tables['D7-5R4La'][PEN+'-5R'])/vol.connectionDistribution['input'][PEN+'-5R']['tot'],\
                                    len(vol.cn_tables['D7-5R4Lb'][PEN+'-5R'])/vol.connectionDistribution['input'][PEN+'-5R']['tot'],\
                                    len(vol.cn_tables['D7-6R3La'][PEN+'-6Ra'])/vol.connectionDistribution['input'][PEN+'-6Ra']['tot'],\
                                    len(vol.cn_tables['D7-6R3La'][PEN+'-6Rb'])/vol.connectionDistribution['input'][PEN+'-6Rb']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb'][PEN+'-6Ra'])/vol.connectionDistribution['input'][PEN+'-6Ra']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb'][PEN+'-6Rb'])/vol.connectionDistribution['input'][PEN+'-6Rb']['tot'],\
                                    len(vol.cn_tables['D7-7R2L'][PEN+'-7R'])/vol.connectionDistribution['input'][PEN+'-7R']['tot'],\
                                    len(vol.cn_tables['D7-2R7L'][PEN+'-7L'])/vol.connectionDistribution['input'][PEN+'-7L']['tot']] ) 
        
        print(mean_PEN_from_D7)
        counts = []
        for D7 in vol.groups['D7']:
            for P in vol.groups[PEN]:
                #discount the ones where D7 outputs in PENs glomerulus
                if not P.neuron_name[5:7] in D7.neuron_name:
                    if vol.connectionDistribution['output'][P.neuron_name]['tot'] > 0:
                        counts.append( len( vol.cn_tables[P.neuron_name][D7.neuron_name] )/vol.connectionDistribution['output'][P.neuron_name]['tot']  )
                    else:
                        print('no connections from ', P.neuron_name)
                        counts.append(0)
                    
        mean_PEN_to_D7 = np.mean(counts) #fractional input FROM PEN
        #print(mean_D7_input_PEN)
        
        for P in vol.groups[PEN]:
            connections_in[P.neuron_name]['D7'] = 5*mean_PEN_from_D7*vol.connectionDistribution['input'][P.neuron_name]['tot'] #~5 D7 arborations per glomerulus
            #the five D7s that input onto P will not receive outputs from P, leaving only 35
            connections_out[P.neuron_name]['D7'] = 35*mean_PEN_to_D7*vol.connectionDistribution['output'][P.neuron_name]['tot'] #assume D7 input uniform
            
        avgCon[PEN]['input']['D7'] = mean_PEN_from_D7 * PEN_in
        #avgCon[PEN]['output']['D7'] = mean_PEN_to_D7 * PEN_out
        
        PENs = [ [PEN+'-6Ra', PEN+'-6Rb'], [PEN+'-5R'], [PEN+'-5L'], [PEN+'-6La', PEN+'-6Lb'] ]
        consout = []
        for i in range(4):
            for EPG in EPGs[i]:
                for P in PENs[i]:
                    consout.append(len(vol.cn_tables[EPG][P]) / vol.connectionDistribution['output'][EPG]['tot'])

        avgCon['EPG']['output'][PEN] = np.mean(consout)*EPG_out
            
        #D7 output is locaized so we take the mean value for pairs of neurons where synapses colocalize
        mean_D7_to_PEN = np.mean( [len(vol.cn_tables['D7-5R4La'][PEN+'-5R'])/vol.connectionDistribution['output']['D7-5R4La']['tot'],\
                                    len(vol.cn_tables['D7-5R4Lb'][PEN+'-5R'])/vol.connectionDistribution['output']['D7-5R4Lb']['tot'],\
                                    len(vol.cn_tables['D7-6R3La'][PEN+'-6Ra'])/vol.connectionDistribution['output']['D7-6R3La']['tot'],\
                                    len(vol.cn_tables['D7-6R3La'][PEN+'-6Rb'])/vol.connectionDistribution['output']['D7-6R3La']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb'][PEN+'-6Ra'])/vol.connectionDistribution['output']['D7-6R3Lb']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb'][PEN+'-6Rb'])/vol.connectionDistribution['output']['D7-6R3Lb']['tot'],\
                                    len(vol.cn_tables['D7-7R2L'][PEN+'-7R'])/vol.connectionDistribution['output']['D7-7R2L']['tot'],\
                                    len(vol.cn_tables['D7-2R7L'][PEN+'-7L'])/vol.connectionDistribution['output']['D7-2R7L']['tot']] ) 
        print(mean_D7_to_PEN)    
    
        counts = []
        for D7 in vol.groups['D7']:
            for P in vol.groups[PEN]:
                if not P.neuron_name[5:7] in D7.neuron_name:
                    #discount P/D7 combinations where D7 outputs in Ps glomerulus; no input in this case
                    counts.append( len( vol.cn_tables[P.neuron_name][D7.neuron_name] )/vol.connectionDistribution['input'][D7.neuron_name]['tot']  )

        mean_D7_from_PEN = np.mean(counts)

        for D7 in vol.groups['D7']:
            connections_out[D7.neuron_name][PEN] = 2*mean_D7_to_PEN*vol.connectionDistribution['output'][D7.neuron_name]['tot'] #each d7 synapse onto 2 PEN
            #d7 does not get input from the two glomeruli where it outputs
            connections_in[D7.neuron_name][PEN] = 14*mean_D7_from_PEN*vol.connectionDistribution['input'][D7.neuron_name]['tot'] #each d7 has 16 PENs synapsing onto it            
            
    mean_EPG_from_D7 = np.mean( [len(vol.cn_tables['D7-5R4La']['EPG-5Ra'])/vol.connectionDistribution['input']['EPG-5Ra']['tot'],\
                                    len(vol.cn_tables['D7-5R4La']['EPG-5Rb'])/vol.connectionDistribution['input']['EPG-5Rb']['tot'],\
                                    len(vol.cn_tables['D7-5R4Lb']['EPG-5Ra'])/vol.connectionDistribution['input']['EPG-5Ra']['tot'],\
                                    len(vol.cn_tables['D7-5R4Lb']['EPG-5Rb'])/vol.connectionDistribution['input']['EPG-5Rb']['tot'],\
                                    len(vol.cn_tables['D7-6R3La']['EPG-6Ra'])/vol.connectionDistribution['input']['EPG-6Ra']['tot'],\
                                    len(vol.cn_tables['D7-6R3La']['EPG-6Rb'])/vol.connectionDistribution['input']['EPG-6Rb']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb']['EPG-6Ra'])/vol.connectionDistribution['input']['EPG-6Ra']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb']['EPG-6Rb'])/vol.connectionDistribution['input']['EPG-6Rb']['tot']] )
    counts = []
    for D7 in vol.groups['D7']:
        for EPG in vol.groups['EPG']:
            if not EPG.neuron_name[4:6] in D7.neuron_name:
                counts.append( len( vol.cn_tables[EPG.neuron_name][D7.neuron_name] )/vol.connectionDistribution['output'][EPG.neuron_name]['tot']  )
    mean_EPG_to_D7 = np.mean(counts)
            
    for EPG in vol.groups['EPG']:
        connections_in[EPG.neuron_name]['D7'] = 5*mean_EPG_from_D7*vol.connectionDistribution['input'][EPG.neuron_name]['tot'] #~4 D7 arborations per glomerulus
        #5 neurons output to EPG, 35 get input
        connections_out[EPG.neuron_name]['D7'] = 35*mean_EPG_to_D7*vol.connectionDistribution['output'][EPG.neuron_name]['tot'] #assume D7 input uniform

    avgCon['EPG']['input']['D7'] = mean_EPG_from_D7 * EPG_in
    avgCon['EPG']['output']['D7'] = mean_EPG_to_D7 * EPG_out 

    mean_D7_to_EPG = np.mean( [len(vol.cn_tables['D7-5R4La']['EPG-5Ra'])/vol.connectionDistribution['input']['D7-5R4La']['tot'],\
                                    len(vol.cn_tables['D7-5R4La']['EPG-5Rb'])/vol.connectionDistribution['output']['D7-5R4La']['tot'],\
                                    len(vol.cn_tables['D7-5R4Lb']['EPG-5Ra'])/vol.connectionDistribution['output']['D7-5R4Lb']['tot'],\
                                    len(vol.cn_tables['D7-5R4Lb']['EPG-5Rb'])/vol.connectionDistribution['output']['D7-5R4Lb']['tot'],\
                                    len(vol.cn_tables['D7-6R3La']['EPG-6Ra'])/vol.connectionDistribution['output']['D7-6R3La']['tot'],\
                                    len(vol.cn_tables['D7-6R3La']['EPG-6Rb'])/vol.connectionDistribution['output']['D7-6R3La']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb']['EPG-6Ra'])/vol.connectionDistribution['output']['D7-6R3Lb']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb']['EPG-6Rb'])/vol.connectionDistribution['output']['D7-6R3Lb']['tot']] )
    counts = []
    for D7 in vol.groups['D7']:
        for EPG in vol.groups['EPG']:
            if not EPG.neuron_name[4:6] in D7.neuron_name:
                counts.append( len( vol.cn_tables[EPG.neuron_name][D7.neuron_name] )/vol.connectionDistribution['input'][D7.neuron_name]['tot']  )
    mean_D7_from_EPG = np.mean(counts)


        
        
    mean_PEG_from_D7 = np.mean( [len(vol.cn_tables['D7-5R4La']['PEG-5R'])/vol.connectionDistribution['input']['PEG-5R']['tot'],\
                                   len(vol.cn_tables['D7-5R4Lb']['PEG-5R'])/vol.connectionDistribution['input']['PEG-5R']['tot'],\
                                    len(vol.cn_tables['D7-6R3La']['PEG-6R'])/vol.connectionDistribution['input']['PEG-6R']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb']['PEG-6R'])/vol.connectionDistribution['input']['PEG-6R']['tot']] )
    counts = []
    for D7 in vol.groups['D7']:
        for PEG in vol.groups['PEG']:
            if not PEG.neuron_name[4:6] in D7.neuron_name:
                if vol.connectionDistribution['output'][PEG.neuron_name]['tot'] > 0:
                    counts.append( len( vol.cn_tables[PEG.neuron_name][D7.neuron_name] ) / vol.connectionDistribution['output'][PEG.neuron_name]['tot']  )
                else:
                    print('no connections from', PEG.neuron_name)
                    counts.append(0)
                #else: counts.append(0)
    mean_PEG_to_D7 = np.mean(counts)
            
    for PEG in vol.groups['PEG']:
        connections_in[PEG.neuron_name]['D7'] = 5*mean_PEG_from_D7*vol.connectionDistribution['input'][PEG.neuron_name]['tot'] #~4 D7 arborations per glomerulus
        connections_out[PEG.neuron_name]['D7'] = 35*mean_PEG_to_D7*vol.connectionDistribution['output'][PEG.neuron_name]['tot'] #assume D7 input uniform

    avgCon['PEG']['input']['D7'] = mean_PEG_from_D7 * PEG_in
    avgCon['PEG']['output']['D7'] = mean_PEG_to_D7 * PEG_out

    mean_D7_to_PEG = np.mean( [len(vol.cn_tables['D7-5R4La']['PEG-5R'])/vol.connectionDistribution['output']['D7-5R4La']['tot'],\
                                   len(vol.cn_tables['D7-5R4Lb']['PEG-5R'])/vol.connectionDistribution['output']['D7-5R4Lb']['tot'],\
                                    len(vol.cn_tables['D7-6R3La']['PEG-6R'])/vol.connectionDistribution['output']['D7-6R3La']['tot'],\
                                    len(vol.cn_tables['D7-6R3Lb']['PEG-6R'])/vol.connectionDistribution['output']['D7-6R3Lb']['tot']] )
    counts = []
    for D7 in vol.groups['D7']:
        for PEG in vol.groups['PEG']:
            if not PEG.neuron_name[4:6] in D7.neuron_name:
                counts.append( len( vol.cn_tables[PEG.neuron_name][D7.neuron_name]/vol.connectionDistribution['input'][D7.neuron_name]['tot'] )  )
    mean_D7_from_PEG = np.mean(counts)
    
        
    countsin = []    
    countsout = [] 
    for D1 in vol.groups['D7']:
        for D2 in vol.groups['D7']:
            countsin.append( len( vol.cn_tables[D1.neuron_name][D2.neuron_name] )/vol.connectionDistribution['input'][D2.neuron_name]['tot']  ) #assume D7s connect uniformly
            countsout.append( len( vol.cn_tables[D1.neuron_name][D2.neuron_name] )/vol.connectionDistribution['output'][D1.neuron_name]['tot']  ) #assume D7s connect uniformly
    mean_D7_in = np.mean(countsin) 
    mean_D7_out = np.mean(countsout) 
        
    for D7 in vol.groups['D7']:
        connections_out[D7.neuron_name]['PEG'] = 2*mean_D7_to_PEG*vol.connectionDistribution['output'][D7.neuron_name]['tot'] #each d7 synapse onto 2 PEG
        #2 PEG get input from D7, 14 give output
        connections_in[D7.neuron_name]['PEG'] = 14*mean_D7_from_PEG*vol.connectionDistribution['input'][D7.neuron_name]['tot'] #each d7 has 16 PEGs synapsing onto it
        connections_out[D7.neuron_name]['EPG'] = 2*56/16*mean_D7_to_EPG*vol.connectionDistribution['output'][D7.neuron_name]['tot'] #each d7 synapse onto 2*56/16 EPG
        #some get input, the rest give output
        connections_in[D7.neuron_name]['EPG'] = (56-56/8)*mean_D7_from_EPG*vol.connectionDistribution['input'][D7.neuron_name]['tot'] #each d7 has 56-56/8 EPGs synapsing onto it
        connections_out[D7.neuron_name]['D7'] = 39*mean_D7_out*vol.connectionDistribution['output'][D7.neuron_name]['tot'] #each d7 synapse onto 39 other D7s
        connections_in[D7.neuron_name]['D7'] = 39*mean_D7_in*vol.connectionDistribution['input'][D7.neuron_name]['tot'] #each d7 has 39 D7s synapsing onto it

   
 
    PEGs = [ ['PEG-6R'], ['PEG-5R'], ['PEG-5L'], ['PEG-6L'] ]
    
    consout = []
    consin = []
    for i in range(4):
        for EPG in EPGs[i]:
            for PEG in PEGs[i]:
                consout.append(len(vol.cn_tables[EPG][PEG]) / vol.connectionDistribution['output'][EPG]['tot'])
                consin.append(len(vol.cn_tables[EPG][PEG]) / vol.connectionDistribution['input'][PEG]['tot'])

    avgCon['EPG']['output']['PEG'] = np.mean(consout)*EPG_out

    avgCon['PEG']['input']['EPG'] = np.mean(consin)*PEG_in
    
    connectionEstimate['input'] = connections_in
    connectionEstimate['output'] = connections_out
    
    return connectionEstimate, avgCon


def get_connection_estimate_NO(vol):
    #we assume that there is zero compartmentalization and we thus only consider left vs right separation
    
    avgCon = {'PEN1':{}, 'PEN2':{}, 'NO-LAL-G':{}}
    for k1, g1 in avgCon.items():
        g1['input'] = {}
        g1['output'] = {}
        
    
    connectionEstimate = {}
        
    connections_in = copy.deepcopy(vol.connectionDistribution['input'])
    connections_out = copy.deepcopy(vol.connectionDistribution['output'])
    
    PEN_avgcounts = {}
    
    for PEN in ['PEN1', 'PEN2']:#we have NO-LAL-Gall in the same noduli as PENx-yR. Consider PEN, PEN2 separately
        
        counts1, counts2 = [], [] #get total inputs/outputs
        for P in vol.groups[PEN]:
            counts1.append(vol.connectionDistribution['input'][P.neuron_name]['tot'])
            counts2.append(vol.connectionDistribution['output'][P.neuron_name]['tot']) 
        PEN_input = np.mean(counts1)
        PEN_output = np.mean(counts2)
        
        PEN_avgcounts[PEN] = [PEN_input, PEN_output] 


    for NO in ['L', 'R']: #now consider left and right nodulus separately for these connection estimates
        
        counts1, counts2 = [], []
        
        #consider interaction of PEN1 with PEN1, PEN2 with PEN2 and PEN1 with PEN2 (both ways)
        for g1 in ['PEN1', 'PEN2']:
            for g2 in ['PEN1', 'PEN2']:
                if g1 == g2: nPartners = 7 #we do not count a self-interactions
                else: nPartners = 8 #each PEN1 interacts with 8 PEN2s and vice versa
                countsin = []
                countsout = []
                for P1 in vol.groups[g1]:
                    for P2 in vol.groups[g2]:
                        if NO in P1.neuron_name and NO in P2.neuron_name: #consider right and left glomerulus individually
                            countsin.append( len( vol.cn_tables[P2.neuron_name][P1.neuron_name] )/vol.connectionDistribution['input'][P1.neuron_name]['tot']  ) #assume PENs connect uniformly
                            countsout.append( len( vol.cn_tables[P1.neuron_name][P2.neuron_name] )/vol.connectionDistribution['output'][P1.neuron_name]['tot']  )
                
                mean_PEN_in = np.mean(countsin) 
                mean_PEN_out = np.mean(countsout)
                print(NO, g1, g2, mean_PEN_in, mean_PEN_out)
                
                avgCon[g1]['input'][g2] = mean_PEN_in*PEN_avgcounts[g1][0]
                avgCon[g1]['output'][g2] = mean_PEN_out*PEN_avgcounts[g1][1]

                for P in vol.groups[g1]:
                    if NO in P.neuron_name:
                        connections_out[P.neuron_name][g2] = nPartners*mean_PEN_out*vol.connectionDistribution['output'][P.neuron_name]['tot'] #each g1 neuron synapses onto npartners g2 neurons
                        connections_in[P.neuron_name][g2] = nPartners*mean_PEN_in*vol.connectionDistribution['input'][P.neuron_name]['tot'] 
                        
        
    connectionEstimate['input'] = connections_in
    connectionEstimate['output'] = connections_out
    
    return connectionEstimate, avgCon

def get_connection_estimate_EB(vol):
    '''
    Need to have a hard think about this...
    
    PEN6R/PEN6L / EPG5L/EPG5R in ta
    PEN5R/PEN7L / EPG6L in tb
    PEN5L/PEN7R / EPG6R in tc
    
    Could consider explicitly... label a, b, c
    
    Divide into the three regions for which we have data.
    Consider interactions within region and to neighboring region. Consider these two types of interactions separately
    Let the ring neurons connect equally to everything: simply do by type over all combinations
    Consider GE to be in wedge 3
    
    How many cross-interactions do we have???
    
    2 GE per wedge, 2 PEN1, 2 PEN2, 2 PEG, 6 EPG
    
    #once I have all this stuff working I should try assuming we have fully traced all PENs, EPGs, PEGs in Wb and see how this compares to the avg number
    
    '''
    
    avgCon = {'GE':{}, 'EPG':{}, 'PEG':{}, 'PEN1':{}, 'PEN2':{}}#, 'R':{}}
    for k1, g1 in avgCon.items():
        g1['input'] = {}
        g1['output'] = {}
    
    connectionEstimate = {}
    connections_in = copy.deepcopy(vol.connectionDistribution['input'])
    connections_out = copy.deepcopy(vol.connectionDistribution['output'])
    
    #Consider first same wedge interactions for EPG
    
    EPGa = ['EPG-6La', 'EPG-6Lb'] #wedge 1
    EPGb = ['EPG-5Ra', 'EPG-5Rb', 'EPG-5Rc', 'EPG-5La', 'EPG-5Lb', 'EPG-5Lc'] #wedge 2
    EPGc = ['EPG-6Ra', 'EPG-6Rb'] #wedge 3
    EPGs = [EPGa, EPGb, EPGc]
    
    PEN1a, PEN1b, PEN1c = ['PEN1-5R', 'PEN1-7L'], ['PEN1-6Ra', 'PEN1-6Rb', 'PEN1-6La', 'PEN1-6Lb'], ['PEN1-7R', 'PEN1-5L']
    PEN1s = [PEN1a, PEN1b, PEN1c]
    PEN2a, PEN2b, PEN2c = ['PEN2-5R', 'PEN2-7L'], ['PEN2-6Ra', 'PEN2-6Rb', 'PEN2-6La', 'PEN2-6Lb'], ['PEN2-7R', 'PEN2-5L']
    PEN2s = [PEN2a, PEN2b, PEN2c]
    
    PEGa, PEGb, PEGc = ['PEG-6L'], ['PEG-5R', 'PEG-5L'], ['PEG-6R']
    PEGs = [PEGa, PEGb, PEGc]
    
    GEa, GEb, GEc = [], [], ['GE-r1']
    GEs = [GEa, GEb, GEc]
    
    Rs = ['Ra', 'Rb', 'Rc', 'Rd', 'Re', 'Rf', 'Rg']
    
    groups = {'EPG': EPGs, 'PEN1': PEN1s, 'PEN2': PEN2s, 'PEG': PEGs, 'GE': GEs}
    
    #first get EPG self interactions in wedge
    
    for label, group in groups.items():
        
        print('considering new group', label) #this is type of neuron
        
        condict = {'to':{'s':{}, 'n':{}}, 'from': {'s':{},'n':{}}} #keep track of connections. s is same, n is neighboring
       
        #first get selfinteractions within wedge
        totin, totout = [], []
        countsin = []
        countsout = []

        for wedge in group: #each 'wedge' is a list of neurons
            for N1 in wedge: #consider all neurons
                totin.append(vol.connectionDistribution['input'][N1]['tot']) #get mean inputs and outputs for this type of neuron
                totout.append(vol.connectionDistribution['output'][N1]['tot'])
                for N2 in wedge: #interacting with all other neurons
                    if N1 != N2: #no self-interaction
                        countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                        countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
                        
        condict['to']['s'][label] = np.mean(countsout) #same wedge
        condict['from']['s'][label] = np.mean(countsin)
        avgCon[label]['input'][label+'_s']=np.mean(countsin)*np.mean(totin)
        avgCon[label]['output'][label+'_s']=np.mean(countsout)*np.mean(totout)
    
        #Now get self-interaction in neighboring wedges
        countsin = []
        countsout = []
        for inds in [ [0,1], [1,0], [1,2], [2,1] ]:
            w1, w2 = group[inds[0]], group[inds[1]] #w1, w2 is now a pair of lists of neurons in neighboring wedges
            for N1 in w1:
                for N2 in w2:
                    countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                    countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
                    
        condict['to']['n'][label] = np.mean(countsout) #neighboring wedge
        condict['from']['n'][label] = np.mean(countsin) #we get a zero for GE-GE since this is empty list, gives nan
        avgCon[label]['input'][label+'_n']=np.mean(countsin)*np.mean(totin)
        avgCon[label]['output'][label+'_n']=np.mean(countsout)*np.mean(totout)
        #now done with all self-interacions. Add them to our connectivity dicts
        
        if label == 'EPG': Ns = 54/8 #number of neurons in a wedge; number of neighboring neurons is Nn = 2Ns
        elif label == 'PEN1' or label == 'PEN2' or label =='PEG' or label =='GE': Ns = 2
        else: print('something wrong, does not recognize', label)
        
        for wedge in group:
            for N in wedge: 
                connections_in[N][label] = \
                    Ns * condict['from']['s'][label] * vol.connectionDistribution['input'][N]['tot'] +\
                    2*Ns * condict['from']['n'][label] * vol.connectionDistribution['input'][N]['tot']
                connections_out[N][label] = \
                    Ns * condict['to']['s'][label] * vol.connectionDistribution['output'][N]['tot'] +\
                    2*Ns * condict['to']['n'][label] * vol.connectionDistribution['output'][N]['tot']
        
        #now do the cross-interactions with neurons from another group
        for label2, group2 in groups.items():
            if label2 != label: #we've already considered self-interactions explicitly
                
                #do same wedge interactions first
                countsin = []
                countsout = []
                for i in range(3):
                    t1s = group[i]
                    t2s = group2[i] #neurons in same wedge of different type
                    
                    for N1 in t1s:
                        for N2 in t2s:
                            countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                            countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
        
                condict['to']['s'][label2] = np.mean(countsout) #same wedge
                condict['from']['s'][label2] = np.mean(countsin)
                avgCon[label]['input'][label2+'_s']=np.mean(countsin)*np.mean(totin)
                avgCon[label]['output'][label2+'_s']=np.mean(countsout)*np.mean(totout)
                
                #now do neighboring wedge interactions
                
                for inds in [ [0,1], [1,0], [1,2], [2,1] ]:
                    t1s = group[inds[0]]
                    t2s = group2[inds[1]]
                    for N1 in t1s:
                        for N2 in t2s:
                            countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                            countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
        
                condict['to']['n'][label2] = np.mean(countsout) #neighboring wedge
                condict['from']['n'][label2] = np.mean(countsin)
                avgCon[label]['input'][label2+'_n']=np.mean(countsin)*np.mean(totin)
                avgCon[label]['output'][label2+'_n']=np.mean(countsout)*np.mean(totout)
                
                #need to add data as appropriate
                #cons_from = 0
                if label2 == 'EPG': Ns = 54/8 #number of neurons in a wedge; number of neighboring neurons is Nn = 2Ns
                elif label2 == 'PEN1' or label2 == 'PEN2' or label2 =='PEG' or label2 =='GE': Ns = 2
                else: print('something wrong, does not recognize', label2)
                
                
                for wedge in group:
                    for N in wedge: 
                        connections_in[N][label2] = \
                                                Ns * condict['from']['s'][label2] * vol.connectionDistribution['input'][N]['tot'] +\
                                                2*Ns * condict['from']['n'][label2] * vol.connectionDistribution['input'][N]['tot']
                        connections_out[N][label2] = \
                                                Ns * condict['to']['s'][label2] * vol.connectionDistribution['output'][N]['tot'] +\
                                                2*Ns * condict['to']['n'][label2] * vol.connectionDistribution['output'][N]['tot']
                
        #finally consider interactions with the rings neurons (we only need input and output from our neurons of interest here). Consider these explicitly
            
        countsin = []
        countsout = []
        for wedge in group:
            for N1 in wedge:
                for R in Rs:
                    countsout.append( len( vol.cn_tables[N1][R] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                    countsin.append( len( vol.cn_tables[R][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
                        
        condict['to']['s']['R'] = np.mean(countsout) #no s/n label necessary but we'll label s for convenience (maintains data structure)
        condict['from']['s']['R'] = np.mean(countsin)
        avgCon[label]['input']['R']=np.mean(countsin)*np.mean(totin)
        avgCon[label]['output']['R']=np.mean(countsout)*np.mean(totout)

        #aaaaand then we add the data to our connection dicts for further use
        print('in', countsout)
        print('out', countsin)
        print('means', 100*np.mean(countsin), 100*np.mean(countsout))
        print('\n')
        for wedge in group:
            for N in wedge: 
                #we assume 100 ring neurons
                connections_in[N]['R'] = 100 * condict['from']['s']['R'] * vol.connectionDistribution['input'][N]['tot']
                connections_out[N]['R'] = 100 * condict['to']['s']['R'] * vol.connectionDistribution['output'][N]['tot']
                
        print(connections_in[N]['R'])
        
    connections_in['GE-r1']['GE']=5*16 #add this manually since it's the only self-interaction. 16 GEs, 5 self-interactions
    connections_out['GE-r1']['GE']=5*16

    connectionEstimate['input'] = connections_in
    connectionEstimate['output'] = connections_out
    
    return connectionEstimate, avgCon


def get_connection_estimate_EB_wedge(vol):
    '''
    
    PEN6R/PEN6L / EPG5L/EPG5R in tilea
    PEN5R/PEN7L / EPG6L in tileb
    PEN5L/PEN7R / EPG6R in tilec
    
    EPG5L ta1, EPG5R ta2
    EPG6L tb1
    EPG6R tc1
    i.e. EPGs in adjacent tiles are not always in adjacent wedges.
    
    Divide into the three tiles for which we have data.
    Consider interactions within tile and to neighboring tile. Consider these two types of interactions separately
    Let the ring neurons connect equally to everything: simply do by type over all combinations
    Consider GE to be in wedge EPG-5R (a2) only adjacent is EPG-5L
    
    How many cross-interactions do we have???
    
    2 GE per tile, 2 PEN1, 2 PEN2, 2 PEG, 6 EPG
    
    Treat everything by tile except for EPG, GE which is treated by wedge...
    
    For GE we will assume that it's in wedge 21, where we don't have any EPG tracing.
    We will assume that the 21:12 GE-EPG ratio is equal to the
    2:1 PEN1 ratio.
    
    In reality this will not be a great assumption since GE overlaps somewhat with w12.
    However, given the lack of tracing and only having 1 GE, it is VERY hard to estimate
    how much of an error we're making and it's not immediately obvious how to do it better
    
    '''
    
    avgCon = {'EPG':{}, 'PEG':{}, 'PEN1':{}, 'PEN2':{}}#, 'R':{}}
    for k1, g1 in avgCon.items():
        g1['input'] = {}
        g1['output'] = {}
    
    connectionEstimate = {}
    connections_in = copy.deepcopy(vol.connectionDistribution['input'])
    connections_out = copy.deepcopy(vol.connectionDistribution['output'])
    
    #Consider first same wedge interactions for EPG
    
    EPGa = ['EPG-6La', 'EPG-6Lb'] #tile 1
    EPGb = ['EPG-5Ra', 'EPG-5Rb', 'EPG-5Rc', 'EPG-5La', 'EPG-5Lb', 'EPG-5Lc'] #tile 2
    EPGc = ['EPG-6Ra', 'EPG-6Rb'] #tile 3
    EPGs = [EPGa, EPGb, EPGc] #at the level of tilles
    
    EPGa1 = ['EPG-6La', 'EPG-6Lb'] #tile 1
    EPGb1 = ['EPG-5La', 'EPG-5Lb', 'EPG-5Lc'] #tile 2
    EPGb2 = ['EPG-5Ra', 'EPG-5Rb', 'EPG-5Rc'] #tile 2
    EPGc2 = ['EPG-6Ra', 'EPG-6Rb'] #tile 3
    EPGws = [EPGa1, EPGb1, EPGb2, EPGc2] #at the level of wedges
    
    
    PEN1a, PEN1b, PEN1c = ['PEN1-5R', 'PEN1-7L'], ['PEN1-6Ra', 'PEN1-6Rb', 'PEN1-6La', 'PEN1-6Lb'], ['PEN1-7R', 'PEN1-5L']
    PEN1s = [PEN1a, PEN1b, PEN1c] #tiles
    PEN2a, PEN2b, PEN2c = ['PEN2-5R', 'PEN2-7L'], ['PEN2-6Ra', 'PEN2-6Rb', 'PEN2-6La', 'PEN2-6Lb'], ['PEN2-7R', 'PEN2-5L']
    PEN2s = [PEN2a, PEN2b, PEN2c] #tiles
    
    PEGa, PEGb, PEGc = ['PEG-6L'], ['PEG-5R', 'PEG-5L'], ['PEG-6R']
    PEGs = [PEGa, PEGb, PEGc] #tiles

    Rs = ['Ra', 'Rb', 'Rc', 'Rd', 'Re', 'Rf', 'Rg']
    
    groups = {'EPG': EPGws, 'PEN1': PEN1s, 'PEN2': PEN2s, 'PEG': PEGs}
    
    #first get EPG self interactions in wedge
    
    totins = {}
    totouts = {}
    
    for label, group in groups.items():
        
        print('considering new group', label) #this is type of neuron
        
        condict = {'to':{'s':{}, 'n':{}}, 'from': {'s':{},'n':{}}} #keep track of connections. s is same, n is neighboring
       
        
        ########################################## Self interactions
        #first get selfinteractions within wedge
        totin, totout = [], []
        countsin = []
        countsout = []

        for wedge in group: #each 'wedge' is a list of neurons
            for N1 in wedge: #consider all neurons
                totin.append(vol.connectionDistribution['input'][N1]['tot']) #get mean inputs and outputs for this type of neuron
                totout.append(vol.connectionDistribution['output'][N1]['tot'])
                
                totins[label] = totin
                totouts[label] = totout
                
                for N2 in wedge: #interacting with all other neurons
                    if N1 != N2: #no self-interaction
                        countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                        countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
                        
        condict['to']['s'][label] = np.mean(countsout) #same wedge
        condict['from']['s'][label] = np.mean(countsin)
        avgCon[label]['input'][label+'_s']=np.mean(countsin)*np.mean(totin)
        avgCon[label]['output'][label+'_s']=np.mean(countsout)*np.mean(totout)
        
        
    
        #Now get self-interaction in neighboring wedges
        countsin = []
        countsout = []
        
        if label == 'EPG': indsN = [ [1,2], [2,1] ] #only tile2 has adjacent wedges
        else: indsN = [ [0,1], [1,0], [1,2], [2,1] ] #consider adjacent tiles for tile neurons
        
        for inds in indsN:
            w1, w2 = group[inds[0]], group[inds[1]] #w1, w2 is now a pair of lists of neurons in neighboring wedges
            for N1 in w1:
                for N2 in w2:
                    countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                    countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
                    
        condict['to']['n'][label] = np.mean(countsout) #neighboring wedge
        condict['from']['n'][label] = np.mean(countsin) #we get a zero for GE-GE since this is empty list, gives nan
        avgCon[label]['input'][label+'_n']=np.mean(countsin)*np.mean(totin)
        avgCon[label]['output'][label+'_n']=np.mean(countsout)*np.mean(totout)
        #now done with all self-interacions. Add them to our connectivity dicts
        
        if label == 'EPG': Ns = 54/16 #number of neurons in a wedge; number of neighboring neurons is Nn = 2Ns
        elif label == 'PEN1' or label == 'PEN2' or label =='PEG' or label =='GE': Ns = 2 #number of neurons in a tile (16/8)
        else: print('something wrong, does not recognize', label)
        
        for wedge in group:
            for N in wedge: 
                connections_in[N][label] = \
                    Ns * condict['from']['s'][label] * vol.connectionDistribution['input'][N]['tot'] +\
                    2*Ns * condict['from']['n'][label] * vol.connectionDistribution['input'][N]['tot']
                connections_out[N][label] = \
                    Ns * condict['to']['s'][label] * vol.connectionDistribution['output'][N]['tot'] +\
                    2*Ns * condict['to']['n'][label] * vol.connectionDistribution['output'][N]['tot']
                    
                    
                    
        ######################################################### cross interactions
        
        #now do the cross-interactions with neurons from another group
        for label2, group2 in groups.items():
            if label2 != label: #we've already considered self-interactions explicitly and EPGs need to be considered by wedge
                
                #do same tile interactions first
                countsin = []
                countsout = []
                
                if label2 == 'EPG': indsS = [ [0,0], [1,1], [1,2], [2,3]  ] #to EPGs same tile
                elif label == 'EPG': indsS = [ [0,0], [1,1], [2,1], [3,2]  ] #from EPGs same tile
                else: indsS = [ [0,0], [1,1], [2,2] ]
                
                
                for inds in indsS:
                    t1s = group[inds[0]] #neurons in tile/wedge of type A
                    t2s = group2[inds[1]] #neurons in same tile/wedge of different type B
                    
                    for N1 in t1s:
                        for N2 in t2s:
                            countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                            countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
        
                condict['to']['s'][label2] = np.mean(countsout) #same wedge
                condict['from']['s'][label2] = np.mean(countsin)
                avgCon[label]['input'][label2+'_s']=np.mean(countsin)*np.mean(totin)
                avgCon[label]['output'][label2+'_s']=np.mean(countsout)*np.mean(totout)
                
                #now do neighboring wedge/tile interactions
                
                if label2 == 'EPG': indsN = [ [0,1], [2,2]  ] #these are the only immediately adjacent wedges
                elif label == 'EPG': indsN = [ [1,0], [2,2] ]
                else: indsN = [ [0,1], [1,0], [1,2], [2,1] ] #pairs of adjacent tiles
                
                for inds in indsN:
                    t1s = group[inds[0]]
                    t2s = group2[inds[1]]
                    for N1 in t1s:
                        for N2 in t2s:
                            countsout.append( len( vol.cn_tables[N1][N2] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                            countsin.append( len( vol.cn_tables[N2][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
        
                condict['to']['n'][label2] = np.mean(countsout) #neighboring wedge/tile
                condict['from']['n'][label2] = np.mean(countsin)
                avgCon[label]['input'][label2+'_n']=np.mean(countsin)*np.mean(totin)
                avgCon[label]['output'][label2+'_n']=np.mean(countsout)*np.mean(totout)
                
                #need to add data as appropriate
                #cons_from = 0
                if label2 == 'EPG':
                    Ns = 54/8 #number of neurons in a tile
                    Nn = Ns #for EPGs we only have Ns neighboring neurons in adjacent wedges
                elif label2 == 'PEN1' or label2 == 'PEN2' or label2 =='PEG' or label2 =='GE':
                    Ns = 2 #number of neurons in a tile
                    Nn = 2*Ns #for tile neurons we ahve twice as many neighboring neurons
                else: print('something wrong, does not recognize', label2)
                
                if label =='EPG': Nn = Nn/2 #EPGs are wedge neurons and only have neighbors on one side
                
                
                
                for wedge in group:
                    for N in wedge: 
                        connections_in[N][label2] = \
                                                Ns * condict['from']['s'][label2] * vol.connectionDistribution['input'][N]['tot'] +\
                                                Nn * condict['from']['n'][label2] * vol.connectionDistribution['input'][N]['tot']
                        connections_out[N][label2] = \
                                                Ns * condict['to']['s'][label2] * vol.connectionDistribution['output'][N]['tot'] +\
                                                Nn * condict['to']['n'][label2] * vol.connectionDistribution['output'][N]['tot']
                
        #finally consider interactions with the rings neurons (we only need input and output from our neurons of interest here). Consider these explicitly
        countsin = []
        countsout = []
        for wedge in group:
            for N1 in wedge:
                for R in Rs:
                    countsout.append( len( vol.cn_tables[N1][R] )/vol.connectionDistribution['output'][N1]['tot']  ) #add output
                    countsin.append( len( vol.cn_tables[R][N1] )/vol.connectionDistribution['input'][N1]['tot']  ) #add input
                        
        condict['to']['s']['R'] = np.mean(countsout) #no s/n label necessary but we'll label s for convenience (maintains data structure)
        condict['from']['s']['R'] = np.mean(countsin)
        avgCon[label]['input']['R']=np.mean(countsin)*np.mean(totin)
        avgCon[label]['output']['R']=np.mean(countsout)*np.mean(totout)

        #aaaaand then we add the data to our connection dicts for further use
        print('in', countsout)
        print('out', countsin)
        print('means', 100*np.mean(countsin), 100*np.mean(countsout))
        print('\n')
        for wedge in group:
            for N in wedge: 
                #we assume 100 ring neurons
                connections_in[N]['R'] = 100 * condict['from']['s']['R'] * vol.connectionDistribution['input'][N]['tot']
                connections_out[N]['R'] = 100 * condict['to']['s']['R'] * vol.connectionDistribution['output'][N]['tot']
                
        print(connections_in[N]['R'])

    
    #avg output to neighboring EPGs is 24.8, input 54.2. average same/neighboring ratio for PEN1, PEN2 is 2.69
            

    connectionEstimate['input'] = connections_in
    connectionEstimate['output'] = connections_out
    
    return connectionEstimate, avgCon
    
    
    


if __name__ == '__main__':
    
    from buildCX import buildCX, get_connection_distributions
    
    Norm = True #determines whether or not to normalize the synapse count
    include_unknown = False #determines whether or not to include an explicit 'unknown' category
    
    show = True #show the neurons as we import them
    vol = buildCX(show=show) #construct a 'volume' instance from the neurons
    
    
    #plot connectivity distribution without extrapolation
    get_connection_distributions(vol, Norm = Norm, include_unknown = include_unknown)
    
    #plot connectivity distributions with extrapolation

    #Protocerebral bridge
    vol.subvolumes['PB'].get_connection_distribution(relation = 'input')
    vol.subvolumes['PB'].get_connection_distribution(relation = 'output') 
    est, avgCon = get_connection_estimate_PB(vol.subvolumes['PB'])
    vol.subvolumes['PB'].connectionEstimate = est
    vol.subvolumes['PB'].plot_connection_distribution(relation = 'input', title='estimate_PB_input', estimate = True, Norm = Norm, include_unknown = include_unknown)
    vol.subvolumes['PB'].plot_connection_distribution(relation = 'output', title='estimate_PB_output', estimate = True, Norm = Norm, include_unknown = include_unknown)
    
    #Noduli
    vol.subvolumes['NO'].get_connection_distribution(relation = 'input')
    vol.subvolumes['NO'].get_connection_distribution(relation = 'output')
    est, avgCon = get_connection_estimate_NO(vol.subvolumes['NO'])
    vol.subvolumes['NO'].connectionEstimate = est
    vol.subvolumes['NO'].plot_connection_distribution(relation = 'input', title='estimate_NO_input', estimate = True, Norm = Norm, include_unknown = include_unknown)
    vol.subvolumes['NO'].plot_connection_distribution(relation = 'output', title='estimate_NO_output', estimate = True, Norm = Norm, include_unknown = include_unknown)
    
    #Ellipsoid body
    vol.subvolumes['EB'].get_connection_distribution(relation = 'input')
    vol.subvolumes['EB'].get_connection_distribution(relation = 'output')
    est, avgCon = get_connection_estimate_EB_wedge(vol.subvolumes['EB'])
    vol.subvolumes['EB'].connectionEstimate = est
    vol.subvolumes['EB'].plot_connection_distribution(relation = 'input', title='estimate_EB_input', estimate = True, Norm = Norm, include_unknown = include_unknown)
    vol.subvolumes['EB'].plot_connection_distribution(relation = 'output', title='estimate_EB_output', estimate = True, Norm = Norm, include_unknown = include_unknown)