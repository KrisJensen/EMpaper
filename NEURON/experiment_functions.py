#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 11:44:07 2018

@author: loaner

functions for NEURON investigations

"""

from neuron import h, gui
import matplotlib.pyplot as plt
import numpy as np
import pickle
import random
plt.rcParams['pdf.fonttype'] = 42

from neuron_simulator import neuron_simulator


def intra_inter_variability(neurons = ['PEN2-5R', 'PEN2-6Ra', 'PEN2-6Rb', 'PEN2-7R',
                                       'PEN2-5L', 'PEN2-6La', 'PEN2-6Lb', 'PEN2-7L'],
                            inputs = ['PEG'], outputs = 'all', regions = ['EB', 'PB', 'NO'],
                            repeat = 10, nsyns = range(1,43), savefig = 'Dan_43syn'):
    '''plots the mean mean and sd of the potential in the different compartments for the neurons in neurons
    injecting inputs at inputs and reading out at outputs.
    PEG only project onto PEN2 in the EB so we can just consider all connections
    when looking at PEG->PEN2 synapses which is convenient.'''
        
    
    if 'PEG' in inputs: Type = 'ACh'
    elif 'D7' in inputs: Type = 'ACh'
    else:
        print 'Type not recognized'
        return
    
    
    allmeans, allsds = [{} for i in range(2)]
    
    for region in regions:
        allmeans[region] = np.zeros((len(neurons), len(nsyns)))
        allsds[region] = np.zeros((len(neurons), len(nsyns)))
        
    for neuroncount, neuron in enumerate(neurons):
        sim = neuron_simulator(h, v_init = -55.0, tstop = 70.0)
        sim.read_neuron('./neurons/'+neuron+'_all.neuron')
        print('  new neuron:', neuron, '  number of inputs:', len(sim.syns_in_tid[inputs[0]]))
        for n, nsyn in enumerate(nsyns):
            newmeans, newsds = [ [[] for i in range(len(regions)) ] for i in range(2) ]
            for rep in range(repeat):
                
                tofire = list(np.random.choice(sim.syns_in_tid[inputs[0]], size=nsyn, replace = False))
                sim.fire(tofire, Type = Type, newstim=[25.0, 3, 2.0, 0.2])#[40.0, 8, 6.0, 0.2])
                #newstim is [interval, number, start, noise]
                
                ns = [0]
                sim.record(['EPG']) #measure of EB output variability
                ns.append(len(sim.v_vecs_record))
                    
                if 'PEG' in inputs:
                    sim.record(sim.syns_in_tid['D7']) #no output in PB so instead measure at D7 input synapses
                elif 'D7' in inputs:
                    sim.record(sim.syns_in_tid['D7']) #measure at EPG inputs
                ns.append(len(sim.v_vecs_record))
        
                if 'NO' in regions:
                    sim.record(['NO-LAL-G']) #measure of NO output variability
                    ns.append(len(sim.v_vecs_record))
                    
                sim.run(Plot = False)#filename = 'PEN2_variability/PEG_EB_PEN2_'+neuron+'_'+region+'_'+str(nsyn))
                    
                for j in range(len(regions)):
                        
                    maxes = []
                    for i, v in enumerate(sim.v_vecs_record[ns[j]:ns[j+1]]): #measure for correct region
                        vec = np.array(v)
                        maxes.append(np.nanmax(vec)) #find max potential for each synapse
                            
                    mean= np.nanmean(maxes)
                    sd = np.nanstd(maxes)
                    newmeans[j].append(mean)
                    newsds[j].append(sd)
                    
                sim.v_vecs_record = [] #reset to allow for new measurements
                sim.i_vec_in = []
                sim.netcons = []
                try:
                    for synapse in tofire: del sim.syns_in[synapse]
                except:
                    continue
                
                print neuron+', # synapses: '+str(nsyn)+', rep:'+str(rep+1)+', mean:'+str(newmeans[0][-1])
                

            for j in range(len(regions)):
                allmeans[regions[j]][neuroncount, n] = np.nanmean(newmeans[j])
                allsds[regions[j]][neuroncount, n] = np.nanmean(newsds[j])
                
    for j in range(len(regions)):
        allmeans[regions[j]] = np.nanmean(allmeans[regions[j]], axis=0)
        allsds[regions[j]] = np.nanmean(allsds[regions[j]], axis=0)
                
        
    pickle.dump([allmeans, allsds], open(savefig+'.pickled', 'wb'))
    
    cols = ['y', 'g', 'b']
    plt.figure(figsize = (5, 3.5))
    for i, region in enumerate(regions):
        plt.plot(nsyns, allmeans[region], '-', color = cols[i])
        plt.fill_between(nsyns, allmeans[region]-allsds[region], allmeans[region]+allsds[region],
                         color=cols[i], alpha=0.25)
    plt.xlabel('Input synapses')
    plt.ylabel('Potential / mV')
    plt.legend(regions)
    plt.savefig(savefig+'.eps', bbox_inches = 'tight')
    plt.savefig(savefig+'.pdf', bbox_inches = 'tight')
    plt.savefig(savefig+'.png', bbox_inches = 'tight', dpi = 240)
    plt.show()
    
    with open(savefig+'.txt', 'w') as f:
        for region in regions:
            f.write(' '.join([str(val) for val in allmeans[region]])+'\n')
            f.write(' '.join([str(val) for val in allsds[region]])+'\n')
    
if __name__ == '__main__':
    
    #plot variability plots
    
    '''
    intra_inter_variability(inputs = ['D7'], repeat = 10, nsyns = range(1,38), savefig = 'results/PEN2_input_PB',
                            neurons = ['PEN2-5R', 'PEN2-6Ra', 'PEN2-6Rb', 'PEN2-7R', 'PEN2-7L'],
                            regions = ['EB', 'PB'])
        
    intra_inter_variability(inputs = ['PEG'], repeat = 10, nsyns = range(1,43), savefig = 'results/PEN2_input_EB',
                            neurons = ['PEN2-5R', 'PEN2-6Ra', 'PEN2-6Rb', 'PEN2-7R', 'PEN2-7L'],
                            regions = ['EB', 'PB'])
    '''
    #plot example plots
    
    
    ##### make variability examples with input in PB
    random.seed(20270606)
    sim = neuron_simulator(h, v_init = -55.0, tstop = 70.0)
    neuron = 'PEN2-7R_all.neuron' #pick example neuron
    sim.read_neuron('./neurons/'+neuron)
    tofire = list(np.random.choice(sim.syns_in_tid['D7-7R2L'], size=30, replace = False)) #pick 30 random synapses
    for i, region in enumerate(['EB', 'PB']):
        sim.fire(tofire, Type = 'ACh', newstim=[25.0, 3, 2.0, 0.2]) #fire an ACh neuron
        if i == 0: sim.record(['EPG'])
        elif i == 1: sim.record(sim.syns_in_tid['D7'])
        sim.run(filename = 'results/example_PEN2_input_PB_record_'+region, write = True, pdf= True)
        sim.v_vecs_record = [] #reset to allow for new measurements
        sim.i_vec_in = []
        sim.netcons = []
        try:
            for synapse in tofire: del sim.syns_in[synapse]
        except:
            continue


        
    

