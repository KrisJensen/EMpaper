#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 15:10:24 2020

@author: kris
"""

import pickle
import copy
from buildCX import buildCX, connectivity, get_connection_distributions
from estimate_synapses_distribution import get_connection_estimate_PB, get_connection_estimate_NO, get_connection_estimate_EB_wedge
from make_PB_EB_extrapolated_matrix import write_combined_cons
from treeNeuron import treeNeuron

show = True #show the neurons as we import them
vol = buildCX(show=show) #construct a 'volume' instance from the neurons


#%% get connectivity matrices
connectivity(vol, show)


#%% get connectivity distributions
Norm = True #determines whether or not to normalize the synapse count
include_unknown = False #determines whether or not to include an explicit 'unknown' category

get_connection_distributions(vol, Norm = Norm, include_unknown = include_unknown)

#%%plot connectivity distributions with extrapolation

#Protocerebral bridge
vol.subvolumes['PB'].get_connection_distribution(relation = 'input')
vol.subvolumes['PB'].get_connection_distribution(relation = 'output') 
est, avgConPB = get_connection_estimate_PB(vol.subvolumes['PB'])
vol.subvolumes['PB'].connectionEstimate = est
vol.subvolumes['PB'].plot_connection_distribution(relation = 'input', title='estimate_PB_input', estimate = True, Norm = Norm, include_unknown = include_unknown)
vol.subvolumes['PB'].plot_connection_distribution(relation = 'output', title='estimate_PB_output', estimate = True, Norm = Norm, include_unknown = include_unknown)

#Noduli
vol.subvolumes['NO'].get_connection_distribution(relation = 'input')
vol.subvolumes['NO'].get_connection_distribution(relation = 'output')
est, avgConNO = get_connection_estimate_NO(vol.subvolumes['NO'])
vol.subvolumes['NO'].connectionEstimate = est
vol.subvolumes['NO'].plot_connection_distribution(relation = 'input', title='estimate_NO_input', estimate = True, Norm = Norm, include_unknown = include_unknown)
vol.subvolumes['NO'].plot_connection_distribution(relation = 'output', title='estimate_NO_output', estimate = True, Norm = Norm, include_unknown = include_unknown)

#Ellipsoid body
vol.subvolumes['EB'].get_connection_distribution(relation = 'input')
vol.subvolumes['EB'].get_connection_distribution(relation = 'output')
est, avgConEB = get_connection_estimate_EB_wedge(vol.subvolumes['EB'])
vol.subvolumes['EB'].connectionEstimate = est
vol.subvolumes['EB'].plot_connection_distribution(relation = 'input', title='estimate_EB_input', estimate = True, Norm = Norm, include_unknown = include_unknown)
vol.subvolumes['EB'].plot_connection_distribution(relation = 'output', title='estimate_EB_output', estimate = True, Norm = Norm, include_unknown = include_unknown)

#store for later use
pickle.dump(avgConEB, open('extrapolated/avg_connections_EB_pickled', 'wb') )
pickle.dump(avgConPB, open('extrapolated/avg_connections_PB_pickled', 'wb') )


#%% create extrapolated connectivity matrices

mod, consEB, consPB = write_combined_cons(continuous = True, fname = 'extrapolated/PB_EB_cont_extrapolated')

#%% parse catmaid neurons for biophysical simulations with NEURON

#### Note this considers all tids independently and takes a while to run ###

mydir = '../NEURON/neurons/'
for neuron in ['PEN2-5R', 'PEN2-6Ra', 'PEN2-6Rb', 'PEN2-7R',
                  'PEN2-5L', 'PEN2-6La', 'PEN2-6Lb', 'PEN2-7L']:
    vol.get_cons_by_neuron()   
    n = copy.copy( vol.get_neuron(neuron) )
    t = treeNeuron(n, vol.neuron_inputs[n.neuron_name], vol.neuron_outputs[n.neuron_name]) #construct a tree
    t.add_synapses(vol) #get synapses from neuropil
    t.write_neuron(mydir+n.neuron_name+'_all') #write to file



