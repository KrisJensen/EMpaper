#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 08:56:24 2018

@author: loaner
"""
import copy
from pymaid import *
import numpy as np
import matplotlib.pyplot as plt
import pickle
plt.rcParams['pdf.fonttype'] = 42


rm = CatmaidInstance('',\
                     '',\
                     '',\
                     '')

class volume():
    '''class for analyzing neurons in the central complex.
    CAn group neurons according to type and split according to neuropil.
    Contains methods for analyzing connectivity'''
    
    def __init__(self, url='', arg1='',\
                 arg2='', token='',\
                 neurons=[], name = 'all'):
        '''initialize by connection to catmaid and setting up data structures.
        Neurons is a list of either skids or names as exact matching strings'''
        
        #rm = CatmaidInstance( url , arg1 , arg2 , token ) #connect to catmaid
        
        self.url, self. arg1, self.arg2, self.token = url, arg1, arg2, token #store for later use
        
        #initialize data structures
        self.groups = dict() #contains subgroups of neurons, e.g. PEN, EPG etc
        self.subvolumes = dict() #contsins subregions e.g. PB, EB, NO
        self.groups['all'] = [] #store neurons as lists for easy iteration        
        self.all_names = dict() #pairs a name with a neuron
        self.all_skids = dict() #pairs a skid with a neuron
        self.neuron_count = 0 #this is also just length of all
        self.cn_tables = {} #connection tables
        self.full_cn_tables_input={} #input tables
        self.full_cn_tables_output={} #output tables
        self.connectionDistribution = {} #distributions
        self.connectionEstimate = {} #extrapolated distributions
        self.neuron_input = {}
        self.neuron_output = {}
        self.all_cons = []
        self.name = name #neuropil name
        
        for neuron in neurons: #initialize with neurons given
            self.add_neuron
        
    def add_neuron(self, neuron, name='default', subgroup = [], to_despike='None', rootnode = 'None', roottype = 'None',\
                   threshold = 'None', show = True):
        '''can take either skid or exact catmaid name as strings, or neuron class isntance.
        Returns neuron. name is the name used for the neuron (string), subgroup is a list of
        any subgroups the neuron should be added to.
        Despiking can help fix broken neurons.'''
        
        if type(neuron) == str or type(neuron) == int:
            neuron = get_neuron(neuron) #fetch from skid or name
            
        if not name == 'default':
            neuron.neuron_name = name #we can rename from the default catmaid names
            
        print('adding neuron', neuron.neuron_name)
            
        if not threshold == 'None':
            despike.despike(neuron, despike, 8000, show = show) # we 'despike' the misaligned nodes for plotting
        elif not to_despike == 'None':
            despike_neuron(neuron, inplace=True, sigma=despike) # this is pymaid's function
            
        self.neuron_count += 1
        self.groups['all'].append(neuron)
        
        neuron.rootnode = rootnode #tid of rootnode for calculating distances
        neuron.roottype = roottype #type of rootnote: 'in' or 'out'

        #maintain dicts of neuron name and skid
        self.all_names[neuron.neuron_name] = neuron #might want to do something if name already exists
        self.all_skids[neuron.skeleton_id] = neuron
        
        if type(subgroup) == list:
            #we're adding to multiple subgroups
            for sg in subgroup:
                if not sg in self.groups.keys(): #initialize if subgroup does not exist
                    self.init_subgroup(sg)   
                self.groups[sg].append(neuron) #add to subgroup
                
        else: #adding to a single subgroup
            if not subgroup in self.groups.keys():
                self.init_subgroup(subgroup)      #initialize if it does not exist      
            self.groups[subgroup].append(neuron)
        
        return neuron
        
    def init_subgroup(self, subgroup, neurons = [ ]):
        '''initiate subgroup of neurons, e.g. corresponding to particular celltype'''
        self.groups[subgroup] = [ ] #add subgroup to dict of subgroups
        for neuron in neurons:
            self.add_to_subgroup(neuron) #add neurons to subgroup
        return self.groups[subgroup]
        
    def add_to_subgroup(self, neuron, subgroup):
        '''add neuron to subgroup'''
        if not subgroup in self.groups:
            self.init_subgroup[subgroup] #initialize if subgroup does not exist
        self.groups[subgroup].append(neuron) #add neurons to subgroup
        
    def subvolume(self, name, neurons=[]):
        '''initiate a new subvolume corresponding to particular neuropil (PB, EB, NO)'''
        skids = [n.skeleton_id for n in neurons] #list of skids to add to subvolume
        #construct new instance of the volume class and add to subvolumes. Add neurons.
        self.subvolumes[name] = volume(self.url, self.arg1, self.arg2, self.token, skids, name = name)

        
    def add_to_subvolume(self, neuron, subvolume, name='default', subgroup=[], rootnode = 'None', roottype = 'None',):
        'add neuron, can be either skid or instance of neuron class'''
        if type(neuron) == str or type(neuron) == int:
            neuron = copy.copy(self.get_neuron(neuron))
        self.subvolumes[subvolume].add_neuron(neuron, name, subgroup, rootnode = rootnode, roottype = roottype)
    
    def split(self, neuron, cut1, cut2, rev=False, reverse=False, show=True, Type ='None'):
        '''split a PEN/EPG neuron into PB, EB and G/NO segments given treenode ids to cut
        'Type' is an ad-hoc parameter that tells us which part of the split corresponds to which neuropil
        this is necessary because we use pymaids cut_neuron function which does not take neuropil identities
        into account'''
        
        neuron = self.get_neuron(neuron)
        print('\ncutting '+neuron.neuron_name)
        
        if Type == 'None': #somewhat ad hoc to cut our specific neurons of interest
            N_PB = cut_neuron(neuron, cut_node=cut1 )[0]
            N_G = cut_neuron(neuron, cut_node=cut2 )[1]
            N_EB = cut_neuron(cut_neuron(neuron, cut_node=cut1 )[1], cut_node=cut2 )[0]
            
        if Type == 'new':
            N_PB = cut_neuron(neuron, cut_node=cut1 )[1]
            N_G = cut_neuron(neuron, cut_node=cut2 )[0]
            N_EB = cut_neuron(cut_neuron(neuron, cut_node=cut1 )[0], cut_node=cut2 )[1]

        if Type == 'rev':
            N_G = cut_neuron(neuron, cut_node=cut1 )[1]
            N_PB = cut_neuron(neuron, cut_node=cut2 )[0]
            N_EB = cut_neuron(cut_neuron(neuron, cut_node=cut1 )[0], cut_node=cut2 )[1]
            
        if Type =='fj': 
            N_PB = cut_neuron(neuron, cut_node=cut1 )[0]
            N_G = cut_neuron(neuron, cut_node=cut2 )[0]
            N_EB = cut_neuron(cut_neuron(neuron, cut_node=cut1 )[1], cut_node=cut2 )[1]

        if Type =='d':
            N_G = cut_neuron(neuron, cut_node=cut2 )[0]
            N_PB = cut_neuron(neuron, cut_node=cut1 )[0]
            N_EB = cut_neuron(cut_neuron(neuron, cut_node=cut1 )[1], cut_node=cut2 )[1]
            
        if reverse: #need to flip the identity
            temp1 = copy.copy(N_G)
            temp2 = copy.copy(N_PB)      
            N_PB = temp1
            N_G = temp2
            
        if show: #plot neuron as cut
            fig, ax = N_PB.plot2d(color='red', method='2d', connectors=False)
            fig, ax = N_EB.plot2d(color='blue', method='2d', connectors=False, ax=ax)
            fig, ax = N_G.plot2d(color='green', method='2d', connectors=False, ax=ax)
            plt.show()

        return N_PB, N_EB, N_G
    
    def get_neuron(self, identifier):
        '''given skid, name or neuron class, returns corresponding instance of neuron class
        allows ua to not worry too much about how the neuron is specified if we run this function before
        subsequent manipulations'''
        if type(identifier) == str:
            try:
                skid = int(identifier)
                neuron = self.all_skids[identifier]
            except:
                neuron = self.all_names[identifier]
        elif type(identifier) == int:
            neuron = self.all_names[ str(identifier) ]
        else:
            try:
                identifier.neuron_name #will raise exception if not neuron type
                neuron = identifier
            except:
                print('identifier not recognized')
                return
        return neuron

    def construct_cn_tables(self, check=False):
        '''constructs dictionary with connections from all neurons in volume to all
        neurons in volume. If check, we check connections against nodes in neurons.
        This is necessary when working with cut neurons.'''
        
        cn_tables = dict() #we make this a dict of dicts of connections cn_tables[pre][post] = connections
        
        for pre in self.groups['all']: #we need distances between all neurons
            pre = self.get_neuron(pre)
            cn_tables[pre.neuron_name] = dict() #make inner dict
            pre_json = self.get_json(pre)
            cn_table = pre_json['outgoing'] #find all connections from the incoming neuron
                
            for post in self.groups['all']:
                
                post = self.get_neuron(post)
                try: connections = cn_table[ post.skeleton_id ]['links'] #get the connections to post
                except KeyError:
                    connections = [] #if there are no connections, the dict is empty
                except ConnectionError:
                    print('No connection, continuing without data') #we might have a network error and not want to crash
                    connections = []

                if check: #if we work with a cut neuron, we remove all entries that aren't in our fragment
                    newconnections = []
                    print(connections)
                    for cn in connections:
                        if cn[0] in pre.nodes.treenode_id.values:
                            print('connection in fragment ', cn[0])
                            newconnections.append(cn) #add to list if in fragment
                    connections = newconnections
                    print(connections)
                
                cn_tables[pre.neuron_name][post.neuron_name]=connections #add to dict
                
        self.cn_tables = cn_tables

    
    def get_neurons(self, group = 'all'):
        return self.groups[group]
    
    def connection_number(self, pre, post, cn_table='None', from_connectors=False):
        '''return number of times pre synapses unto post. Pre, post are neurons
        can take either skids, names or neurons. if from_connectors is true, we do this from a list of connectors
        This works for neuron fragments which from_partners does not. Could also use our cn_tables but wrote this first'''
        
        pre = self.get_neuron(pre)
        skid = self.get_neuron(post).skeleton_id
        
        if type(cn_table) == str: #if we already have a cn table we can skip this
            if from_connectors: cn_table = cn_table_from_connectors(pre) #get connection table for fragment
            elif cn_table == 'stored': return len(self.cn_tables[pre][post])
            else: cn_table = pre.partners #get connection table from intact neuron

        downstream = cn_table[ cn_table.relation=='downstream' ] #we want to look at what this synapses onto

        try: nSyn = downstream[ downstream.skeleton_id == skid ].total.values[0] #this has number of synapses
        except IndexError: nSyn=0 #if there are no synapses
        except TypeError: 
            try: nSyn = downstream[ downstream.skeleton_id == int(skid) ].total.values[0] #inconsistencies in pymaid means we sometimes need to use an int as key
            except IndexError: nSyn = 0

        return nSyn        
    
    def allconnections(self, group1='all', group2='all', file='None', method='w', table = 'none',\
                       heatmap='None', from_connectors = True, show=True, grouping='default', size = (15,15)):
        '''given two groups of neurons, give matrix of all connections FROM group1 TO group2
        file specifies name of file to write data to, method whether we want to create new file or append it
        heatmap gives name of file to write heatmap to. From connectors specifies whether or not we're using a frgment
        or an intact neuron'''
        
        if group1 == 'all': group1 = self.groups['all']
        if group2 == 'all': group2 = self.groups['all']
        data=[]
        tot=0.0
        
        if grouping == 'EB':
            g1 = ['PEN1-5R', 'PEN1-7L', 'PEN2-5R', 'PEN2-7L', 'EPG-6La', 'EPG-6Lb', 'PEG-6L',\
                      'PEN1-6Ra', 'PEN1-6Rb', 'PEN1-6La', 'PEN1-6Lb','PEN2-6Ra', 'PEN2-6Rb', 'PEN2-6La', 'PEN2-6Lb',\
                      'EPG-5Ra', 'EPG-5Rb', 'EPG-5Rc', 'EPG-5La', 'EPG-5Lb', 'EPG-5Lc', 'PEG-5R', 'PEG-5L',\
                      'PEN1-7R', 'PEN1-5L', 'PEN2-7R', 'PEN2-5L', 'EPG-6Ra', 'EPG-6Rb', 'PEG-6R',\
                      'Ra','Rb','Rc','Rd','Re','Rf','Rg']
            group1 = [self.get_neuron(n) for n in g1]
            group2 = copy.copy(group1)
        
        for pre in group1: #consider each presynaptic neuron      
            if not table == 'stored':
                if from_connectors: cn_table = cn_table_from_connectors(pre)
                else: cn_table = pre.partners                 
            data.append([])
            for post in group2:
                if table == 'stored':
                    data[-1].append(len(self.cn_tables[pre.neuron_name][post.neuron_name]))
                else:
                    data[-1].append( self.connection_number(pre, post, cn_table=cn_table) ) #find connections to post
            tot += np.sum(data[-1]) #also keep track of total connectivity
        if not file == 'None': #we write our result to a file
            with open('connectivity_matrices/'+file+'.txt', method) as f:
                f.write('{:<10}'.format(' '))
                for post in group2: f.write( '{:<10}'.format( post.neuron_name[:10] )+' ' )
                f.write( '\n' )
                for i, pre in enumerate(group1):
                    f.write( '{:<10}'.format( pre.neuron_name[:10] )+' ' )
                    for j, post in enumerate(group2):
                        f.write( '{:<10}'.format( str(data[i][j])[:10] )+' ' )
                    f.write('\n')
                    
        if not heatmap == 'None': #make a heatmap
            xlabel  = [ post.neuron_name for post in group2 ]
            ylabel = [ pre.neuron_name for pre in group1 ]
            xticks = len(xlabel)
            yticks = len(ylabel)
            
            self.heatmap(data, xticks=xticks, yticks=yticks, xlabel=xlabel, ylabel=ylabel,\
                         title=heatmap, show=show, size = size)
            
        return tot, data
    
    def get_total_connections(self, neuron, relation, Type = 'local'):
        '''return the total number of input connections to this neuron'''
        
        if relation == 'input': relation = 'upstream'
        if relation == 'output': relation = 'downstream'
        
        neuron = self.get_neuron(neuron)
        if Type == 'global':
            cn_table = cn_table_from_connectors(neuron)
            pickle.dump(cn_table, open('pickled_neurons/neuron_connections/'+self.name+'_'+relation+'_'+neuron.neuron_name+'_pickled', 'wb'))
        else:
            cn_table = pickle.load(open('pickled_neurons/neuron_connections/'+self.name+'_'+relation+'_'+neuron.neuron_name+'_pickled', 'rb'))
        connections = cn_table[ cn_table.relation==relation ] #get all INPUTS/OUTPUTS for a neuron
        
        if relation == 'upstream': self.full_cn_tables_input[neuron.neuron_name] = connections
        if relation == 'downstream': self.full_cn_tables_output[neuron.neuron_name] = connections 
        nconnections = sum(connections.total.values)
        
        return nconnections
        

    def get_connection_distribution(self, relation = 'input'):
        '''find the distribution of inputs to each type of neuron from the other types of neurons
        surely this doesn't really work since we don't have the full connectivity of all of our neurons?'''
        
        connections = {}
        for neuron in self.groups['all']:
            print('new neuron', neuron.neuron_name)
            name = neuron.neuron_name
            connections[name] = {}
            connections[name]['tot'] = self.get_total_connections(neuron, relation)
            for group in self.groups.keys():
                if group != 'all':
                    cons = 0
                    for partner in self.groups[group]:
                        if relation == 'input': cons += len(self.cn_tables[partner.neuron_name][name])
                        if relation == 'output': cons += len(self.cn_tables[name][partner.neuron_name])
                    connections[name][group]=cons
            print(connections[name])
        self.connectionDistribution[relation] = connections
            
        return connections
        
    
    def plot_connection_distribution(self, relation = 'input', neurons_init='all',\
                                     title='connection_distribution', show=True, estimate=True,
                                     Norm = True, include_unknown = False):
        '''
        plot distribution of input/output connections as given by 'relation'
        if estimate, we do this for extrapolated synapse counts
        '''

        if neurons_init == 'all': neurons_init = self.groups['all']
        neurons = []
        #first we normalize the connection distribution per neuron
        if estimate: connections = copy.deepcopy(self.connectionEstimate[relation])
        else: connections = copy.deepcopy(self.connectionDistribution[relation])
    
        groups = ['PEN1', 'PEN2', 'EPG', 'PEG', 'D7']
        leg = ['PEN1 (8)', 'PEN2 (8)', 'EPG (10)', 'PEG (4)', 'D7 (7)']
        if estimate: leg = ['PEN1 (16)', 'PEN2 (16)', 'EPG (54)', 'PEG (16)', 'D7 (40)']
        colors = [ '#f1c40f', '#cd6155','#27ae60', '#2e86c1', '#c2c2c2']
                  
        if 'R' in self.groups.keys():
            groups.append('R')
            if estimate: leg.append('R (100)')
            else: leg.append('R (7)')
            colors.append( '#d7bde2' )
                          
        if include_unknown: #include a category for unknown partners
            colors.append('#e8e8e8')
            leg.append('unknown')
            groups.append('unknown')
                          
        for i, neuron in enumerate(neurons_init):
            if int(connections[neuron.neuron_name]['tot']) < 10:
                print('fewer than 10 connections for '+neuron.neuron_name+', not included in plot')

            else:
                if not (neuron.neuron_name in ['Ra','Rb','Rc','Rd','Re','Rf','Rg']):
                    cum = 0.0
                    for key, number in connections[neuron.neuron_name].items():
                        #print(number, self.connectionDistribution[neuron.neuron_name]['tot'])
                        if Norm: norm = number / connections[neuron.neuron_name]['tot']
                        else: norm = number
                        
                        #print(norm)
                        if not key =='tot':
                            connections[neuron.neuron_name][key] = norm
                            if key in groups: cum += norm
                    neurons.append(neuron)
                    if include_unknown:
                        if Norm: connections[neuron.neuron_name]['unknown'] = 1-cum
                        else: connections[neuron.neuron_name]['unknown'] = connections[neuron.neuron_name]['tot']-cum
        N = len(neurons)
        
        fname = 'connectivity_distributions/'+title
        if include_unknown: fname += '_u'
        if not Norm: fname += '_raw'
        f = open(fname+'.txt', 'w') #write to a file as well
        vals = []
        for group in groups:
            if group != 'all':
                vals.append([])
                for neuron in neurons:
                    vals[-1].append( connections[neuron.neuron_name][group] )
                    f.write(group+' '+neuron.neuron_name+' '+str(connections[neuron.neuron_name][group])+'\n')   
        f.close()
        ind = np.arange(N)    # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence
        ps = []
        bottom = np.array( [ 0.0 for x in range(N) ] )
        
        #plot the results
        fig, ax = plt.subplots(figsize=(10,10))
        for i, val in enumerate(vals):
            print('bottom: ', bottom)
            val = np.maximum(val, 0)
            ps.append( plt.bar(ind, val, width, color=colors[i], bottom=bottom) )
            bottom += np.array(val)
        print(vals)
        plt.title(title)
        ax.set_xticks(ind)
        labels = [ ( neuron.neuron_name+' ('+str(self.connectionDistribution[relation][neuron.neuron_name]['tot'])+')' )\
                  for neuron in neurons ]
        ax.set_xticklabels( labels )
        plt.setp(ax.get_xticklabels(), rotation=-90, ha="right", va='top')
        if Norm: plt.yticks(np.arange(0, 1.1, 0.1))
        
        plt.legend([ p[0] for p in ps ],  leg )

        plt.savefig(fname+'.pdf', dpi=360, bbox_inches = 'tight', transparent=True)

        if show: plt.show()
    
        
    def heatmap(self, data,  xticks=None, yticks=None, xlabel=None,\
                ylabel=None, title=None, vmin=None, vmax=None, show='False', size = (15,15) ):
        '''construct heatmap of connectivity matrices'''
        
        fig, ax = plt.subplots(figsize=(15,15))
        im = ax.imshow(data, cmap='coolwarm', vmin=vmin, vmax=vmax) #generate heatmap; could add more options
        cb = fig.colorbar(im, shrink=1.1, aspect=5 ) #add colorbar
        cb.ax.tick_params(labelsize = size[0]/1.5)

        ax.set_xticks(np.arange(xticks))
        ax.set_yticks(np.arange(yticks))
        ax.set_xticklabels(xlabel) #add cellnames
        ax.set_yticklabels(ylabel) #add cellnames
        
        ax.xaxis.set_label_position('top') 
        ax.xaxis.tick_top()
        # Rotate the tick labels and set their alignment
        plt.setp(ax.get_xticklabels(), rotation=-45, ha="right", va='bottom', rotation_mode="anchor") #make it look nice

        plt.xlabel(title, fontsize = int(size[0]*1.5))
        plt.tight_layout()
        print('saving heatmap')
        plt.savefig('connectivity_matrices/'+title+'.png', dpi=360)
        if show: plt.show()
        
    def plot2D(self, neurons='all', color = False, show = True, save='None'):
        '''plots a neuron in 2 dimensions using the pymaid plot2d function'''
        
        if neurons == 'all': neurons = self.groups['all']
        else: neurons = [self.get_neuron(n) for n in neurons] #neurons to plot
        if type(color) == list: color=color
        elif type(color) == str: color = [color for x in range(len(neurons))]
        else: color = ['red' for i in range(len(neurons))] #can specify colors of neurons
        for i, neuron in enumerate(neurons): #plot all of the neurons on the same set of axes
            try: fig, ax = neuron.plot2d(color=color[i], method='2d',\
                                         connectors=False, ax=ax)
            except:fig, ax = neuron.plot2d(color=color[i], method='2d', connectors=False)
        if not save == 'None': #save the figure
            plt.savefig(save, dpi=360, bbox_inches = 'tight')
        if show: plt.show()

    def get_cons_by_neuron(self):
        inp = {}
        out = {}
        tot = []
        for neuron in self.groups['all']:
            inp[neuron.neuron_name] = []
            out[neuron.neuron_name] = []
        for n1, dic in self.cn_tables.items():
            for n2, cons in dic.items():
                for con in cons:
                    out[n1].append( con[0] )
                    inp[n2].append( con[1] )
                    tot.append( ( con[0], con[1] ) )
        self.neuron_inputs = inp
        self.neuron_outputs = out
        self.all_cons = tot

