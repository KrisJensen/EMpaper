#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 21:51:00 2018

@author: loaner

numerical simulation

Gouwens & Wilson j. neuroscience 2009
E_rest = -60 mV
Rm = 20800 Ohm cm2
Cm = 0.79 uF cm-2
Ra (Ri) = 266 ohm cm

http://www.jneurosci.org/content/19/13/5311 (drosophila embryos, 3-9 divisions)
Ionotropic cholinergic mEPSPs
reversal potential near 0 mV (embryonic neurons)
G = I/V = 16/75 = 16/75000 uS (amplitude -16 pA at -75 mV)
#rise time 0.61 ms
#decay time 2.11 ms

http://www.jneurosci.org/content/19/13/5311 (drosophila embryos, 3-9 divisions)
Ionotropic cholinergic mEPSPs
reversal potential near 0 mV (embryonic neurons)
G = I/V = 17.5/85 = 17.5/85000 uS (amplitude -16 pA at -75 mV)
#rise time 0.61 ms
#decay time 2.11 ms

https://www.physiology.org/doi/full/10.1152/jn.2002.88.2.847
Ionotropic GABAergic receptors (Rdl)
G = I/V = -100 pA/(-70-(-55)) = 100/15000 uS
rise time 1 ms
decay time 50 ms
reversal potential -55 mV

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1175349/
GABA-B channels: 6-8 ms lag G-protein activity
gmax = 50/14000
tau1 = 100 ms
tau2 = 195 ms
https://www.frontiersin.org/files/Articles/152407/fncom-09-00139-HTML/image_m/fncom-09-00139-t001.jpg
Er = Ek = -72

"""
#from mpi4py import MPI
from neuron import h, gui
import matplotlib.pyplot as plt
#import neuron
import numpy as np
import time


class neuron_simulator():
    
    '''class for parsing neurons based on EM data (from .neuron files written by treeNeuron.py)
    and seeting up NEURON simulations, inputting according to fitted synaptic inputs and reading outputs
    at specified nodes. This class deals with single neurons only.'''
    
    def __init__(self, h=h, v_init = -55.0, tstop = 100.0):
        '''specify length of simulation (tstop, ms) and resting potential (v_init)'''
        
        h.v_init = v_init #set initial potential to resting potential
        self.h = h
        self.h.tstop = tstop
        
        self.stims = {} #stores NetStim objects
        self.v_vecs_record = []
        self.stimlist = []
        self.gmax = {'ACh':17.5/85000.0, 'GABA':100.0/15000.0, 'GABA-B':50.0/14000.0, 'AP':16000.0/75000.0 } #max conductance of different types of synapses
        self.delays = {'ACh':0.0, 'GABA':0.0, 'GABA-B':10.0, 'AP':0.0 } #10 ms delay for GABA-B from G protein action
        self.cols = []
        
    def reset(self):
        '''reset inputs and recordings for new experiment. Otherwise hoc will still include vectors etc.'''
        self.syns_in, self.syns_out, self.i_vecs_in, self.v_vecs_out, self.netcons, self.v_vecs_record, self.stimlist =\
            {}, {}, [], [], [], [], []
        
    
    def read_neuron(self, filename, fire = [], fire_ids = 'all', record = [], record_ids = 'all', Print = False):
        '''Given a .neuron file, reads in the geometry of the neuron to constuct a NEURON object. Also stores locations
        of all synapses and rootnodes. If 'fire' and 'record' are given (as a neuron name or type), these
        will be set up for firing/recording
        
        each line in the .neuron file is: section# length diameter parent_section# treenode_id
        roots are given by: %root xxx yyy zzz
        inputs by: %post xxx yyy zzz www
        outputs by: %pre xxx yyy zzzz www
        
        comments must have # as first character in the line
        
        '''
        
        #default neurons to fire
        if not type(fire) == list: fire = [fire]
        if not type(record) == list: record = [record]
        
        syns_in_tid = {} #tids of incoming synapses
        syns_in = {} #construct dicts of synapses labelled by celltype
        syns_out_tid = {}
        syns_in_pos = {} #gives us the segments at which we have synapses for later recording
        syns_out = {} #dict of outgoing synapses
        syns_out_pos = {} #positions of outgoing synapses
        i_vecs_in = [] #store current through input synapses
        v_vecs_out = [] #store potential at output synapses
        sections = [] #store sections
        treenode_ids = [] #list of treenode ids by section number
        treenode_dict = {} #dict of section number by treenode_id
        syns_in_pos_dict = {}
        syns_out_pos_dict = {}
        
        #can play around with putting this elsewhere, although order of make netstim, make syn, make con IS important
        stim = self.h.NetStim()
        stim.interval = 20.0 #time between spikes
        stim.number = 4 #number of spikes
        stim.start = 5.0 
        stim.noise = 0.0 #interval is interval*(1-noise)+exp_distribution(u = noise)
        
        self.stims['default'] = stim #can stimulate multiple synapses with this
        self.v_vecs_record = []
        
        with open(filename, 'r') as f:
            for line in f:
                if not (line[0] == '#' or line[0] == '%'): # '#' are comments, '%' special lines
                    split = line.split()
                    if Print: print(split)
                    treenode_dict[split[4]] = len(treenode_ids) #zero indexed
                    treenode_ids.append(split[4])
                    sections.append(self.h.Section(name = split[4])) #name according to tid
                    sections[-1].L = float(split[1]) #length in microns
                    sections[-1].diam = float(split[2]) #diameter in microns
                    sections[-1].nseg=int( max(1, float(split[1])/2.0))#25.0 )) #can get away with must larger segments, but most of them are smaller anyways
                    if not split[3] == '-1': #root has parent -1
                        sections[-1].connect( sections[ int(split[3]) ], 1, 0) #attach 0-end to parent's 1-end (0 indexed)
                        
                if line[0] == '%':
                    if line[0:5] == '%root': # '%root xxx yyy zzz'
                        roots = []
                        split = line.split()
                        syns_out_pos['root'] = []
                        self.roots = []
                        root_str = ''
                        for root in split[1:len(split)]:#iterate over tids
                            roots.append(treenode_dict[root])#store tid
                            syns_out_pos['root'].append( treenode_dict[root] )#store position
                            root_str = root_str+str(treenode_dict[root]) + '_'
                            self.roots.append(int(root))
                        
                    elif line[0:5] == '%post': # '%post PEG 17 19 50'
                        split = line.split()
                        syns_in[split[1]] = []
                        syns_in_pos[split[1]] = []
                        syns_in_tid[split[1]] = []
                        if len(split) > 2:

                            for synid in split[2:]:
                                sec = treenode_dict[synid] # section#
                                syns_in_pos[split[1]].append(sec) #store position
                                syns_in_pos_dict[int(synid)] = sec #call section from treenode_id
                                syns_in_tid[split[1]].append(int(synid)) #store treenodeid vs input neuron name/type
                                
                    elif line[0:4] == '%pre': #'%pre EPG 7 45 56'
                        split = line.split()
                        syns_out[split[1]] = []
                        syns_out_pos[split[1]] = []
                        syns_out_tid[split[1]] = []
                        if len(split) > 2:
                            for synid in split[2:]:
                                #make synapse and set parameters
                                sec = treenode_dict[synid]
                                syns_out_pos[split[1]].append(sec)#position of synapse
                                syns_out_pos_dict[int(synid)] = sec
                                syns_out_tid[split[1]].append(int(synid))

                    elif line[0:7] == '%neuron': #can store the name of the neuron we're looking at
                        self.neuron = line.split()[1]
                                    
        
        #set passive membrane parameters according to results in Gouwens & WIlson
        self.h('forall Ra = 266') #set axial resistivity in ohm cm
        self.h('forall cm = 0.79') #set specific capacitance in uF/cm2
        self.h('forall insert pas') #make a model using passive conductance
        self.h('forall g_pas = %.10f' % (1./20800.) ) #g_pas (conductance) in S/cm2 - 1 over Rm in ohm cm2
        self.h('forall e_pas = %.10f' % (self.h.v_init) ) #e_pas in mV          

        self.sections = sections #need to save these for NEURON to work
        self.syns_in_tid=syns_in_tid
        self.syns_out_tid=syns_out_tid
            
        netcons = []
        for group in fire:
            print(syns_in)
            print 'firing '+group+' #neuron: '+str(len(syns_in[group]))
            if fire_ids == 'all': #can subset which neurons we want to fire. This isn't that useful after implementation of the 'fire' function
                syns = syns_in[group]
                syns_pos = syns_in_pos[group]
            else:
                syns = [syns_in[group][i] for i in fire_ids]
                syns_pos = [syns_in_pos[group][i] for i in fire_ids]
            for pos in syns_pos:
                syn = self.get_syn(pos) #construct ACh synapse as default, can make this flexible later
                syns_in[group].append(syn) #store synapse
                i_vec = self.h.Vector()
                i_vec.record(syn._ref_i)
                i_vecs_in.append(i_vec) #record current through synapse
                netcons.append(self.h.NetCon(stim, syn)) #construct connection to allow us to stimulate this
                netcons[-1].weight[0]=self.gmax['ACh'] #store netcon
                
        for group in record:
            if record_ids == 'all': syns = syns_out_pos[group]
            else: syns = [syns_out_pos[group][i] for i in record_ids]
            for pos in syns:
                v_vec = self.h.Vector()
                v_vec.record(sections[pos](1.0)._ref_v)
                v_vecs_out.append(v_vec) #record potential at these positions
            
        v_vecs_root = []    
        for root in roots:    
            v_vec_root = self.h.Vector()
            v_vec_root.record(sections[root](1.0)._ref_v)
            v_vecs_root.append(v_vec_root) #record potential at roots

        self.t_vec = self.h.Vector()
        self.t_vec.record(self.h._ref_t) #record time

        self.sections, self.syns_in, self.syns_out, self.i_vecs_in, self.v_vecs_out, self.v_vecs_root,\
            self.netcons, self.stim, self.treenode_dict,\
            self.syns_in_pos_dict, self.syns_out_pos_dict, self.syns_in_pos, self.syns_out_pos = \
            sections, syns_in, syns_out, i_vecs_in, v_vecs_out, v_vecs_root, netcons, stim,\
            treenode_dict, syns_in_pos_dict, syns_out_pos_dict, syns_in_pos, syns_out_pos
            
    
    def fire(self, nodes, stim = 'default', Type = 'ACh', gmax = 'default', newstim ='None', sources = [], thresholds = [], segments = []):
        '''give a list of notes that should fire in a subsequent record. List as treenode_ids
        if newstim is not none, we create a new stimulus with [interval, number, start, noise]'''
    
        if not type(nodes) == list: nodes = [nodes]
        if not type(sources) == list: sources = [sources]
        if not type(thresholds) == list: thresholds = [thresholds]
        if not type(segments) == list: segments = [segments]
    
        if len(sources) == 0:
            stim = self.stims[stim] #specifies which stimulus to use. In time, we should allow this to be other neurons
    
        for i, n in enumerate(nodes): #for each synapse
            
            if not newstim == 'None': #given paramters for new stimuli
                stim = self.h.NetStim()
                stim.interval, stim.number, stim.start, stim.noise = newstim
                self.stimlist.append(stim)
                
                
            syn = self.get_syn(self.treenode_dict[str(n)], Type=Type) #construct synapse of Type
            self.syns_in[n] = syn
            
            i_vec = self.h.Vector() #record current through synapse
            i_vec.record(syn._ref_i)
            self.i_vecs_in.append(i_vec)
            
            if len(sources) > 0:
                print(sources[i])
                print(syn)
                h('access '+str(sources[i]))  
                self.netcons.append(self.h.NetCon(sources[i](1.0)._ref_v, syn))
                self.netcons[-1].threshold = thresholds[i]
            
            else: self.netcons.append(self.h.NetCon(stim, syn)) #connect synapse to NetStim
            
            if gmax == 'default': self.netcons[-1].weight[0]=self.gmax[Type] #specified by type
            else: self.netcons[-1].weight[0]=gmax #we can artificially alter gmax if we have a reason to do so
            
            self.netcons[-1].delay = self.delays[Type] #metabotropic receptors often have some inherent delay
    
    def record(self, nodes, segment=1.0, col = 'None'):
        '''nodes must be a list of either INTS of treenode ids, or STRINGS of groups of groups of neurons;
        i.e. [223532, 46456, 2342] or ['EPG', 'PEG', 'PEN1-5R']'''
        
        if not type(nodes) == list: nodes = [nodes]
        
        count = 0 #count number of recorded synapses
        
        for n in nodes:
            
            if type(n) == str:
                for pos in self.syns_out_pos[n]: #we can give group of recording synapses, e.g. 'EPG-5R'
                    v_vec = self.h.Vector()
                    v_vec.record(self.sections[pos](segment)._ref_v) #record potential
                    self.v_vecs_record.append(v_vec)
                    count += 1
                
            else:
                v_vec = self.h.Vector()
                v_vec.record(self.sections[self.treenode_dict[str(n)]](segment)._ref_v) #record potential
                self.v_vecs_record.append(v_vec)
                count += 1
                
        if not col == 'None': self.cols += [col for x in range(count)] #can specify a color in which to plot this group of positions
    
    
    def show(self, t_vec, v_vecs, filename = 'None', title = 'None', legend = 'None', ylim = 'None', col = 'None',\
             ylabel = 'potential / mV', xlabel = 'time / ms', write = False, pdf=False):
        '''plots all of v_vecs (list) vs. time on the same plots. Can specify a list of colors for plotting by col'''
        
        
        plt.figure(figsize = (8,8))
        for i in range(len(v_vecs)):
            
            if len(v_vecs[i]) > 0: #only plot vectors with content
                if col == 'None': plt.plot(t_vec, v_vecs[i], ':')
                else: plt.plot(t_vec, v_vecs[i], ':', color = col[i]) #can optionally specify a color
        
        if not ylim == 'None': plt.ylim(ylim) #can specify ylim
            
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        if not title == 'None': plt.title(title)
        
        if not legend == 'None':
            plt.legend(legend)
        
        if not filename == 'None': #we save fig if a filename is specified
            plt.savefig(filename+'.png', dpi=360)
            if pdf: plt.savefig(filename+'.pdf')
            
            if write:
                with open(filename+'.txt', 'w') as f:
                    for i in range(len(v_vecs)):
                        f.write(' '.join([str(val) for val in t_vec])+'\n')
                        if len(v_vecs[i]) > 0: #only plot vectors with content
                            f.write(' '.join([str(val) for val in np.array(v_vecs[i])])+'\n')
        plt.show()
        
    def run(self, filename = 'None', title = 'None', col = 'None', legend = 'None', Plot = True,
            write = False, pdf=False):
        '''runs a simulation given parameters of h and plots results
        (input current, root potentials and vec_record potentials)'''
        
        t0 = time.time()
        self.h.run() #this takes a while
        dt = time.time()-t0
        print('\n\nsimulation took '+str(np.round(dt, 2))+' seconds')
        
        if Plot:
            self.show(self.t_vec, self.i_vecs_in, title = 'Inputs', filename = filename+'_inputs',\
                      ylabel = 'synaptic current / nA', write=write, pdf=pdf)
            
            
            self.show(self.t_vec, self.v_vecs_record, title = 'Outputs', filename = filename+'_outputs',\
                      ylim = [-72,0], col = col, legend = legend, write=write, pdf=pdf)
            
            
            self.show(self.t_vec, self.v_vecs_root, title = 'Root',\
                      filename = filename+'_root', ylim = [-72,0], write=write, pdf=pdf)#, legend = ['PB', 'EB', 'EB', 'NO'])

    
    def runsim(self, filename = 'PEN2-5R.neuron', fire = [], record = [], fire_ids = 'all', record_ids = 'all', imgfile = 'None'):
        '''let's quickly run a simulation for a neuron specified by filename, firing fire and recording record'''

        self.read_neuron(filename, fire = fire, record = record, fire_ids = fire_ids, record_ids=record_ids)
        self.h.run()
        self.show(self.t_vec, self.v_vecs_in, title = 'Inputs', filename = imgfile+'_inputs')
        self.show(self.t_vec, self.v_vecs_out, title = 'Outputs', filename = imgfile+'_outputs')
        self.show(self.t_vec, self.v_vecs_root, title = 'Root', filename = imgfile+'_root', legend = ['PB', 'EB', 'NO'])
        
        for vec in self.v_vecs_root:
            print(max(list(vec))) #print the max root potentials

    def get_syn(self, pos, Type = 'ACh'):
        '''given a synapse type, returns an Exp2Syn instance with parameters to roughly match experimental data
        for this synapse'''
        syn = self.h.Exp2Syn(self.sections[pos](1.0))
        if Type == 'ACh':
            syn.tau1 = 0.61
            syn.tau2 = 2.11
            syn.e = 0.0
        elif Type == 'GABA':
            syn.tau1 = 5.00
            syn.tau2 = 48.0
            syn.e = -55.0
        elif Type == 'GABA-B':
            syn.tau1 = 60.0
            syn.tau2 = 140.0
            syn.e = -72.0
        elif Type == 'AP': #rapid peak
            syn.tau1 = 0.01
            syn.tau2 = 0.30
            syn.e = 0.0
        return syn
    
