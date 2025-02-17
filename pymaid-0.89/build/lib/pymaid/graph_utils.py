#    This script is part of pymaid (http://www.github.com/schlegelp/pymaid).
#    Copyright (C) 2017 Philipp Schlegel
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along

""" Collection of tools to maniuplate CATMAID neurons using Graph
representations.
"""

import pandas as pd
import numpy as np
import networkx as nx

from scipy.sparse import csgraph, csr_matrix

from pymaid import graph, core, utils, config

# Set up logging
logger = config.logger

__all__ = sorted(['classify_nodes', 'cut_neuron', 'longest_neurite',
                  'split_into_fragments', 'reroot_neuron', 'distal_to',
                  'dist_between', 'find_main_branchpoint',
                  'generate_list_of_childs', 'geodesic_matrix',
                  'subset_neuron', 'node_label_sorting'])


def _generate_segments(x, weight=None):
    """ Generate segments maximizing segment lengths.

    Parameters
    ----------
    x :         CatmaidNeuron | CatmaidNeuronList
                May contain multiple neurons.
    weight :    'weight' | None, optional
                If ``"weight"`` use physical length to determine segment
                length. If ``None`` use number of nodes.

    Returns
    -------
    list
                Segments as list of lists containing treenode ids. List is
                sorted by segment lengths.
    """

    if isinstance(x, pd.DataFrame) or isinstance(x, core.CatmaidNeuronList):
        return [_generate_segments(x.loc[i],
                                   weight=weight) for i in range(x.shape[0])]
    elif isinstance(x, core.CatmaidNeuron):
        pass
    else:
        logger.error('Unexpected datatype: %s' % str(type(x)))
        raise ValueError

    if weight == 'weight':
        # Get distances from end nodes to root
        m = geodesic_matrix(x,
                            directed=True,
                            weight=weight,
                            tn_ids=x.nodes[x.nodes.type == 'end'].treenode_id.values)

        # Sort by distance to root
        endNodeIDs = m.sort_values(x.root[0], ascending=False).index.values
    elif not weight:
        d = _edge_count_to_root(x)
        endNodeIDs = x.nodes[x.nodes.type == 'end'].treenode_id.values
        endNodeIDs = sorted(endNodeIDs, key=d.get, reverse=True)
    else:
        raise ValueError('Unable to use weight "{}"'.format(weight))

    if x.igraph and config.use_igraph:
        g = x.igraph
        # Convert endNodeIDs to indices
        id2ix = {n: ix for ix, n in zip(g.vs.indices,
                                        g.vs.get_attribute_values('node_id'))}
        endNodeIDs = [id2ix[n] for n in endNodeIDs]
    else:
        g = x.graph

    seen = set()
    sequences = []
    for nodeID in endNodeIDs:
        sequence = [nodeID]
        parents = list(g.successors(nodeID))
        while True:
            if not parents:
                break
            parentID = parents[0]
            sequence.append(parentID)
            if parentID in seen:
                break
            seen.add(parentID)
            parents = list(g.successors(parentID))

        if len(sequence) > 1:
            sequences.append(sequence)

    # If igraph, turn indices back to node IDs
    if x.igraph and config.use_igraph:
        ix2id = {v: k for k, v in id2ix.items()}
        sequences = [[ix2id[ix] for ix in s] for s in sequences]

    # Sort sequences by length
    if weight == 'weight':
        sequences = sorted(
            sequences, key=lambda x: m.loc[x[0], x[-1]], reverse=True)
    else:
        sequences = sorted(
            sequences, key=lambda x: d[x[0]] - d[x[-1]], reverse=True)

    return sequences


def _break_segments(x):
    """ Break neuron into small segments connecting ends, branches and root.

    Parameters
    ----------
    x :         CatmaidNeuron | CatmaidNeuronList
                May contain multiple neurons.

    Returns
    -------
    list
                Segments as list of lists containing treenode ids.

    """

    if isinstance(x, pd.DataFrame) or isinstance(x, core.CatmaidNeuronList):
        return [_break_segments(x.loc[i]) for i in range(x.shape[0])]
    elif isinstance(x, core.CatmaidNeuron):
        pass
    else:
        logger.error('Unexpected datatype: %s' % str(type(x)))
        raise ValueError

    if x.igraph and config.use_igraph:
        g = x.igraph
        end = g.vs.select(_indegree=0).indices
        branch = g.vs.select(_indegree_gt=1, _outdegree=1).indices
        root = g.vs.select(_outdegree=0).indices
        seeds = branch + end
        # Converting to set speeds up the "parent in stops" check
        stops = set(branch + root)
        seg_list = []
        for s in seeds:
            parent = g.successors(s)[0]
            seg = [s, parent]
            while parent not in stops:
                parent = g.successors(parent)[0]
                seg.append(parent)
            seg_list.append(seg)
        # Translate indices to treenode IDs
        ix_id = {v: n for v, n in zip(g.vs.indices,
                                      g.vs.get_attribute_values('node_id'))}
        seg_list = [[ix_id[n] for n in s] for s in seg_list]
    else:
        seeds = x.nodes[x.nodes.type.isin(['branch', 'end'])].treenode_id.values
        stops = x.nodes[x.nodes.type.isin(['branch', 'root'])].treenode_id.values
        # Converting to set speeds up the "parent in stops" check
        stops = set(stops)
        g = x.graph
        seg_list = []
        for s in seeds:
            parent = next(g.successors(s), None)
            seg = [s, parent]
            while parent not in stops:
                parent = next(g.successors(parent), None)
                seg.append(parent)
            seg_list.append(seg)

    return seg_list


def _edge_count_to_root(x):
    """ Return a map of nodeID vs number of edges from the first node that
    lacks successors (aka the root).
    """
    if x.igraph and config.use_igraph:
        g = x.igraph
        current_level = g.vs(_outdegree=0).indices
    else:
        g = x.graph
        current_level = list(x.root)

    dist = {}
    count = 1
    next_level = []
    while current_level:
        # Consume all elements in current_level
        while current_level:
            node = current_level.pop()
            dist[node] = count
            next_level.extend(g.predecessors(node))
        # Rotate lists (current_level is now empty)
        current_level, next_level = next_level, current_level
        count += 1

    # Map vertex index to treenode ID
    if x.igraph and config.use_igraph:
        dist = {x.igraph.vs[k]['node_id']: v for k, v in dist.items()}

    return dist


def classify_nodes(x, inplace=True):
    """ Classifies neuron's treenodes into end nodes, branches, slabs
    or root. Adds ``'type'`` column to ``x.nodes``.

    Parameters
    ----------
    x :         CatmaidNeuron | CatmaidNeuronList
                Neuron(s) whose nodes to classify nodes.
    inplace :   bool, optional
                If ``False``, nodes will be classified on a copy which is then
                returned leaving the original neuron unchanged.

    Returns
    -------
    CatmaidNeuron/List
                Copy of original neuron. Only if ``inplace=False``.

    """

    if not inplace:
        x = x.copy()

    # If more than one neuron
    if isinstance(x, (pd.DataFrame, core.CatmaidNeuronList)):
        for i in config.trange(x.shape[0], desc='Classifying'):
            classify_nodes(x.ix[i], inplace=True)
    elif isinstance(x, (pd.Series, core.CatmaidNeuron)):
        # Make sure there are nodes to classify
        if x.nodes.shape[0] != 0:
            if x.igraph and config.use_igraph:
                # Get graph representation of neuron
                vs = x.igraph.vs
                # Get branch/end nodes based on their degree of connectivity
                ends = vs.select(_indegree=0).get_attribute_values('node_id')
                branches = vs.select(
                    _indegree_gt=1).get_attribute_values('node_id')
            else:
                # Get graph representation of neuron
                g = x.graph
                # Get branch/end nodes based on their degree of connectivity
                deg = pd.DataFrame.from_dict(dict(g.degree()), orient='index')
                # [ n for n in g.nodes if g.degree(n) == 1 ]
                ends = deg[deg.iloc[:, 0] == 1].index.values
                # [ n for n in g.nodes if g.degree(n) > 2 ]
                branches = deg[deg.iloc[:, 0] > 2].index.values

            x.nodes.loc[:, 'type'] = 'slab'
            x.nodes.loc[x.nodes.treenode_id.isin(ends), 'type'] = 'end'
            x.nodes.loc[x.nodes.treenode_id.isin(branches), 'type'] = 'branch'
            x.nodes.loc[x.nodes.parent_id.isnull(), 'type'] = 'root'
    else:
        raise TypeError('Unknown neuron type: %s' % str(type(x)))

    if not inplace:
        return x


def distal_to(x, a=None, b=None):
    """ Checks if nodes A are distal to nodes B.

    Important
    ---------
    Please note that if node A is not distal to node B, this does **not**
    automatically mean it is proximal instead: if nodes are on different
    branches, they are neither distal nor proximal to one another! To test
    for this case run a->b and b->a - if both return ``False``, nodes are on
    different branches.

    Also: if a and b are the same node, this function will return ``True``!

    Parameters
    ----------
    x :     CatmaidNeuron
    a,b :   single treenode ID | list of treenode IDs | None, optional
            If no treenode IDs are provided, will consider all treenodes.

    Returns
    -------
    bool
            If *a* and *b* are single treenode IDs respectively.
    pd.DataFrame
            If *a* and/or *b* are lists of treenode IDs. Columns and rows
            (index) represent treenode IDs. Neurons *a* are rows, neurons *b*
            are columns.

    Examples
    --------
    >>> # Get a neuron
    >>> x = pymaid.get_neuron(16)
    >>> # Get a treenode ID from tag
    >>> a = x.tags['TEST_TAG'][0]
    >>> # Check all nodes if they are distal or proximal to that tag
    >>> df = pymaid.distal_to(x, a)
    >>> # Get the IDs of the nodes that are distal
    >>> df[ df[a] ].index.values

    """

    if isinstance(x, core.CatmaidNeuronList) and len(x) == 1:
        x = x[0]

    if not isinstance(x, core.CatmaidNeuron):
        raise ValueError('Please pass a SINGLE CatmaidNeuron')

    if not isinstance(a, type(None)):
        a = utils._make_iterable(a)
        # Make sure we're dealing with integers
        a = np.unique(a).astype(int)
    else:
        a = x.nodes.treenode_id.values

    if not isinstance(b, type(None)):
        b = utils._make_iterable(b)
        # Make sure we're dealing with integers
        b = np.unique(b).astype(int)
    else:
        b = x.nodes.treenode_id.values

    if x.igraph and config.use_igraph:
        # Map treenodeID to index
        id2ix = {n: v for v, n in zip(x.igraph.vs.indices,
                                      x.igraph.vs['node_id'])}

        # Convert treenode IDs to indices
        if isinstance(a, type(None)):
            a = x.igraph.vs.indices
        else:
            a = [id2ix[n] for n in a]
        if isinstance(b, type(None)):
            b = x.igraph.vs
        else:
            b = [id2ix[n] for n in b]

        # Get path lengths
        le = x.igraph.shortest_paths(a, b, mode='OUT')

        # Converting to numpy array first is ~2X as fast
        le = np.asarray(le)

        # Convert to True/False
        le = le != float('inf')

        df = pd.DataFrame(le,
                          index=x.igraph.vs[a]['node_id'],
                          columns=x.igraph.vs[b]['node_id'])
    else:
        # Generate empty DataFrame
        df = pd.DataFrame(np.zeros((len(a), len(b)), dtype=bool),
                          columns=b, index=a)

        # Iterate over all targets
        for nB in config.tqdm(b, desc='Querying paths',
                       disable=(len(b) < 1000) | config.pbar_hide,
                       leave=config.pbar_leave):
            # Get all paths TO this target. This function returns a dictionary:
            # { source1 : path_length, source2 : path_length, ... } containing
            # all nodes distal to this node.
            paths = nx.shortest_path_length(x.graph, source=None, target=nB)
            # Check if sources are among our targets
            df[nB] = [nA in paths for nA in a]

    if df.shape == (1, 1):
        return df.values[0][0]
    else:
        # Return boolean
        return df


def geodesic_matrix(x, tn_ids=None, directed=False, weight='weight'):
    """ Generates geodesic ("along-the-arbor") distance matrix for treenodes
    of given neuron.

    Parameters
    ----------
    x :         CatmaidNeuron | CatmaidNeuronList
                If list, must contain a SINGLE neuron.
    tn_ids :    list | numpy.ndarray, optional
                Treenode IDs. If provided, will compute distances only FROM
                this subset to all other nodes.
    directed :  bool, optional
                If True, pairs without a child->parent path will be returned
                with ``distance = "inf"``.
    weight :    'weight' | None, optional
                If ``weight`` distances are given as physical length.
                if ``None`` distances is number of nodes.

    Returns
    -------
    pd.SparseDataFrame
                Geodesic distance matrix. Distances in nanometres.

    See Also
    --------
    :func:`~pymaid.distal_to`
        Check if a node A is distal to node B.
    :func:`~pymaid.dist_between`
        Get point-to-point geodesic distances.
    """

    if isinstance(x, core.CatmaidNeuronList):
        if len(x) == 1:
            x = x[0]
        else:
            raise ValueError('Cannot process more than a single neuron.')
    elif isinstance(x, core.CatmaidNeuron):
        pass
    else:
        raise ValueError(
            'Unable to process data of type "{0}"'.format(type(x)))

    if x.igraph and config.use_igraph:
        nodeList = x.igraph.vs.get_attribute_values('node_id')

        # Matrix is ordered by vertex number
        m = _igraph_to_sparse(x.igraph, weight_attr=weight)
    else:
        nodeList = tuple(x.graph.nodes())

        m = nx.to_scipy_sparse_matrix(x.graph, nodeList,
                                      weight=weight)

    if not isinstance(tn_ids, type(None)):
        tn_ids = set(utils._make_iterable(tn_ids))
        tn_indices = tuple(i for i, node in enumerate(
            nodeList) if node in tn_ids)
        ix = [nodeList[i] for i in tn_indices]
    else:
        tn_indices = None
        ix = nodeList

    dmat = csgraph.dijkstra(m,
                            directed=directed, indices=tn_indices)

    return pd.SparseDataFrame(dmat, columns=nodeList, index=ix,
                              default_fill_value=float('inf'))


def dist_between(x, a, b):
    """ Returns the geodesic distance between treenodes in nanometers.

    Parameters
    ----------
    x :             CatmaidNeuron | CatmaidNeuronList
                    Neuron containing the nodes.
    a,b :           treenode IDs
                    Treenodes to check.

    Returns
    -------
    int
                    distance in nm

    See Also
    --------
    :func:`~pymaid.distal_to`
        Check if a node A is distal to node B.
    :func:`~pymaid.geodesic_matrix`
        Get all-by-all geodesic distance matrix.

    """

    if isinstance(x, core.CatmaidNeuronList):
        if len(x) == 1:
            x = x[0]
        else:
            raise ValueError('Need a single CatmaidNeuron, got {}'.format(len(x)))

    if isinstance(x, core.CatmaidNeuron):
        if x.igraph and config.use_igraph:
            g = x.igraph
        else:
            g = x.graph
    elif isinstance(x, nx.DiGraph):
        g = x
    elif 'igraph' in str(type(x.igraph)):
        # We can't use isinstance here because igraph library might not be installed
        g = x
    else:
        raise ValueError('Unable to process data of type {0}'.format(type(x)))

    if (utils._is_iterable(a) and len(a) > 1) or \
       (utils._is_iterable(b) and len(b) > 1):
        raise ValueError('Can only process single treenodes. Use pymaid.geodesic_matrix instead.')

    a = utils._make_non_iterable(a)
    b = utils._make_non_iterable(b)

    try:
        _ = int(a)
        _ = int(b)
    except:
        raise ValueError('a, b need to be treenode IDs!')

    # If we're working with network X DiGraph
    if isinstance(g, nx.DiGraph):
        return int(nx.algorithms.shortest_path_length(g.to_undirected(as_view=True),
                                                      a, b,
                                                      weight='weight'))
    else:
        # If not, we're assuming g is an iGraph object
        return g.shortest_paths(g.vs.find(node_id=a),
                                g.vs.find(node_id=b),
                                weights='weight',
                                mode='ALL')[0][0]


def find_main_branchpoint(x, reroot_to_soma=False):
    """ Returns the branch point at which the two largest branches converge.

    Parameters
    ----------
    x :                 CatmaidNeuron | CatmaidNeuronList
                        May contain multiple neurons.
    reroot_to_soma :    bool, optional
                        If True, neuron will be rerooted to soma.

    Returns
    -------
    treenode ID

    """

    # Make a copy
    x = x.copy()

    if isinstance(x, core.CatmaidNeuronList) and len(x) > 1:
        return np.array([find_main_branchpoint(n, reroot_to_soma=reroot_to_soma) for n in x])
    elif isinstance(x, core.CatmaidNeuronList) and len(x) == 1:
        x = x[0]
    elif not isinstance(x, (core.CatmaidNeuron, core.CatmaidNeuronList)):
        raise TypeError(
            'Must provide CatmaidNeuron/List, not "{0}"'.format(type(x)))

    g = graph.neuron2nx(x)

    # First, find longest path
    longest = nx.dag_longest_path(g)

    # Remove longest path
    g.remove_nodes_from(longest)

    # Find second longst path
    sc_longest = nx.dag_longest_path(g)

    # Parent of the last node in sc_longest is the common branch point
    bp = list(x.graph.successors(sc_longest[-1]))[0]

    return bp


def split_into_fragments(x, n=2, min_size=None, reroot_to_soma=False):
    """ Splits neuron into fragments.

    Notes
    -----
    Cuts are based on longest neurites: the first cut is made where the second
    largest neurite merges onto the largest neurite, the second cut is made
    where the third largest neurite merges into either of the first fragments
    and so on.

    Parameters
    ----------
    x :                 CatmaidNeuron | CatmaidNeuronList
                        May contain only a single neuron.
    n :                 int, optional
                        Number of fragments to split into. Must be >1.
    min_size :          int, optional
                        Minimum size of fragment in um to be cut off. If too
                        small, will stop cutting. This takes only the longest
                        path in each fragment into account!
    reroot_to_soma :    bool, optional
                        If True, neuron will be rerooted to soma.

    Returns
    -------
    CatmaidNeuronList

    Examples
    --------
    >>> x = pymaid.get_neuron('16')
    >>> # Cut into two fragments
    >>> cut1 = pymaid.split_into_fragments(x, n=2)
    >>> # Cut into fragments of >10 um size
    >>> cut2 = pymaid.split_into_fragments(x, n=float('inf'), min_size=10)

    """

    if isinstance(x, core.CatmaidNeuron):
        pass
    elif isinstance(x, core.CatmaidNeuronList):
        if x.shape[0] == 1:
            x = x[0]
        else:
            logger.error(
                '%i neurons provided. Please provide only a single neuron!' % x.shape[0])
            raise Exception
    else:
        raise TypeError('Unable to process data of type "{0}"'.format(type(x)))

    if n < 2:
        raise ValueError('Number of fragments must be at least 2.')

    if reroot_to_soma and x.soma:
        x.reroot(x.soma)

    # Collect treenodes of the n longest neurites
    tn_to_preserve = []
    fragments = []
    i = 0
    while i < n:
        if tn_to_preserve:
            # Generate fresh graph
            g = graph.neuron2nx(x)

            # Remove nodes that we have already preserved
            g.remove_nodes_from(tn_to_preserve)
        else:
            g = x.graph

        # Get path
        longest_path = nx.dag_longest_path(g)

        # Check if fragment is still long enough
        if min_size:
            this_length = sum([v / 1000 for k, v in nx.get_edge_attributes(
                x.graph, 'weight').items() if k[1] in longest_path])
            if this_length <= min_size:
                break

        tn_to_preserve += longest_path
        fragments.append(longest_path)

        i += 1

    # Next, make some virtual cuts and get the complement of treenodes for each fragment
    graphs = [x.graph.copy()]
    for fr in fragments[1:]:
        this_g = nx.bfs_tree(x.graph, fr[-1], reverse=True)

        graphs.append(this_g)

    # Next, we need to remove treenodes that are in subsequent graphs from those graphs
    for i, g in enumerate(graphs):
        for g2 in graphs[i + 1:]:
            g.remove_nodes_from(g2.nodes)

    # Now make neurons
    nl = core.CatmaidNeuronList(
        [subset_neuron(x, g, clear_temp=True) for g in graphs])

    # Rename neurons
    for i, n in enumerate(nl):
        n.neuron_name += '_{}'.format(i)

    return nl


def longest_neurite(x, n=1, reroot_to_soma=False, inplace=False):
    """ Returns a neuron consisting of only the longest neurite(s) based on
    geodesic distance.

    Parameters
    ----------
    x :                 CatmaidNeuron | CatmaidNeuronList
                        May contain only a single neuron.
    n :                 int | slice, optional
                        Number of longest neurites to preserve. For example:
                         - ``n=1`` preserves the longest neurites
                         - ``n=2`` preserves two longest neurites
                         - ``n=slice(1,None)`` removes the longest neurite
    reroot_to_soma :    bool, optional
                        If True, neuron will be rerooted to soma.
    inplace :           bool, optional
                        If False, copy of the neuron will be trimmed down to
                        longest neurite and returned.

    Returns
    -------
    CatmaidNeuron
                        Pruned neuron. Only if ``inplace==False``.

    See Also
    --------
    :func:`~pymaid.split_into_fragments`
            Split neuron into fragments based on longest neurites.

    """

    if isinstance(x, core.CatmaidNeuron):
        pass
    elif isinstance(x, core.CatmaidNeuronList):
        if x.shape[0] == 1:
            x = x[0]
        else:
            raise ValueError('Please provide only a single neuron.')
    else:
        raise TypeError('Unable to process data of type "{0}"'.format(type(x)))

    if n < 1:
        raise ValueError(
            'Number of longest neurites to preserve must be at least 1.')

    if not inplace:
        x = x.copy()

    if reroot_to_soma and x.soma:
        x.reroot(x.soma)

    segments = _generate_segments(x, weight='weight')

    if isinstance(n, (int, np.int)):
        tn_to_preserve = [tn for s in segments[:n] for tn in s]
    elif isinstance(n, slice):
        tn_to_preserve = [tn for s in segments[n] for tn in s]
    else:
        raise TypeError('Unable to use N of type "{}"'.format(type(n)))

    subset_neuron(x, tn_to_preserve, inplace=True)

    if not inplace:
        return x


def reroot_neuron(x, new_root, inplace=False):
    """ Reroot neuron to new root.

    Parameters
    ----------
    x :        CatmaidNeuron | CatmaidNeuronList
               List must contain a SINGLE neuron.
    new_root : int | str
               Node ID or tag of the node to reroot to.
    inplace :  bool, optional
               If True the input neuron will be rerooted.

    Returns
    -------
    CatmaidNeuron
               Rerooted neuron. Only if ``inplace=False``.

    See Also
    --------
    :func:`~pymaid.CatmaidNeuron.reroot`
                Quick access to reroot directly from CatmaidNeuron/List
                objects.

    Examples
    --------
    >>> n = pymaid.get_neuron(16)
    >>> # Reroot neuron to its soma
    >>> n2 = pymaid.reroot_neuron(n, n.soma)

    """

    if new_root is None:
        raise ValueError('New root can not be <None>')

    if isinstance(x, core.CatmaidNeuron):
        pass
    elif isinstance(x, core.CatmaidNeuronList):
        if x.shape[0] == 1:
            x = x.loc[0]
        else:
            raise Exception(
                '{0} neurons provided. Please provide only a single neuron!'.format(x.shape[0]))
    else:
        raise Exception('Unable to process data of type "{0}"'.format(type(x)))

    # If new root is a tag, rather than a ID, try finding that node
    if isinstance(new_root, str):
        if new_root not in x.tags:
            logger.error(
                '#{}: Found no treenodes with tag {} - please double check!'.format(x.skeleton_id, new_root))
            return
        elif len(x.tags[new_root]) > 1:
            logger.error(
                '#{}: Found multiple treenodes with tag {} - please double check!'.format(x.skeleton_id, new_root))
            return
        else:
            new_root = x.tags[new_root][0]

    if not inplace:
        x = x.copy()

    # Skip if new root is old root
    if x.root == new_root:
        if not inplace:
            return x
        else:
            return

    if x.igraph and config.use_igraph:
        path = x.igraph.get_shortest_paths(x.igraph.vs.find(node_id=new_root),
                                           x.igraph.vs.find(node_id=x.root))[0]
        epath = x.igraph.get_shortest_paths(x.igraph.vs.find(node_id=new_root),
                                            x.igraph.vs.find(node_id=x.root),
                                            output='epath')[0]

        edges = [(s, t) for s, t in zip(path[:-1], path[1:])]

        weights = [x.igraph.es[e]['weight'] for e in epath]

        # Get all weights and append inversed new weights
        all_weights = x.igraph.es['weight'] + weights

        # Add inverse edges: old_root->new_root
        x.igraph.add_edges([(e[1], e[0]) for e in edges])

        # Re-set weights
        x.igraph.es['weight'] = all_weights

        # Remove new_root->old_root
        x.igraph.delete_edges(edges)

        # Get degree of old root for later categorisation
        old_root_deg = len(x.igraph.es.select(_target=path[-1]))

        # Translate path indices to treenode IDs
        ix2id = {ix: n for ix, n in zip(x.igraph.vs.indices,
                                        x.igraph.vs.get_attribute_values('node_id'))}
        path = [ix2id[i] for i in path]
    else:
        # If this NetworkX graph is just an (immutable) view, turn it into a
        # full, independent graph
        if isinstance(x.graph, nx.classes.graphviews.ReadOnlyGraph):
            x.graph = nx.DiGraph(x.graph)
        g = x.graph

        # Walk from new root to old root and remove edges along the way
        parent = next(g.successors(new_root), None)
        if not parent:
            # new_root is already the root
            return
        path = [new_root]
        weights = []
        while parent is not None:
            weights.append(g[path[-1]][parent]['weight'])
            g.remove_edge(path[-1], parent)
            path.append(parent)
            parent = next(g.successors(parent), None)

        # Invert path and add weights
        new_edges = [(path[i + 1], path[i],
                      {'weight': weights[i]}) for i in range(len(path) - 1)]

        # Add inverted path between old and new root
        g.add_edges_from(new_edges)

        # Get degree of old root for later categorisation
        old_root_deg = g.in_degree(path[-1])

    # Propagate changes in graph back to treenode table
    x.nodes.set_index('treenode_id', inplace=True)
    x.nodes.loc[path[1:], 'parent_id'] = path[:-1]
    if old_root_deg == 1:
        x.nodes.loc[path[-1], 'type'] = 'slab'
    elif old_root_deg > 1:
        x.nodes.loc[path[-1], 'type'] = 'branch'
    else:
        x.nodes.loc[path[-1], 'type'] = 'end'
    x.nodes.reset_index(drop=False, inplace=True)

    # Set new root's parent to None
    x.nodes.parent_id = x.nodes.parent_id.astype(object)
    x.nodes.loc[x.nodes.treenode_id == new_root, 'parent_id'] = None

    if x.igraph and config.use_igraph:
        x._clear_temp_attr(exclude=['igraph', 'classify_nodes'])
    else:
        x._clear_temp_attr(exclude=['graph', 'classify_nodes'])

    if not inplace:
        return x
    else:
        return


def cut_neuron(x, cut_node, ret='both'):
    """ Split neuron at given point and returns two new neurons.

    Note
    ----
    Split is performed between cut node and its parent node. However, cut node
    will still be present in both resulting neurons.

    Parameters
    ----------
    x :        CatmaidNeuron | CatmaidNeuronList
               Must be a single neuron.
    cut_node : int | str | list
               Node ID(s) or a tag(s) of the node(s) to cut. Multiple cuts are
               performed in the order of ``cut_node``. Fragments are ordered
               distal -> proximal.
    ret :      'proximal' | 'distal' | 'both', optional
               Define which parts of the neuron to return. Use this to speed
               up processing when making only a single cut!

    Returns
    -------
    distal -> proximal :    CatmaidNeuronList
                            Distal and proximal part of the neuron. Only if
                            ``ret='both'``. The distal->proximal order of
                            fragments is tried to be maintained for multiple
                            cuts but is not guaranteed.
    distal :                CatmaidNeuronList
                            Distal part of the neuron. Only if
                            ``ret='distal'``.
    proximal :              CatmaidNeuronList
                            Proximal part of the neuron. Only if
                            ``ret='proximal'``.

    Examples
    --------
    First example: single cut at node tag

    >>> import pymaid
    >>> rm = pymaid.CatmaidInstance(url, http_user, http_pw, token)
    >>> n = pymaid.get_neuron(16)
    >>> # Cut neuron
    >>> nl = cut_neuron(n, 'SCHLEGEL_LH')
    >>> nl
    <class 'pymaid.core.CatmaidNeuronList'> of 2 neurons
                     neuron_name skeleton_id  n_nodes  n_connectors  \
    0  PN glomerulus VA6 017 DB          16     3603          1295
    1  PN glomerulus VA6 017 DB          16     9142           741
       n_branch_nodes  n_end_nodes  open_ends  cable_length review_status   soma
    0             303          323          3    960.118516            NA  False
    1             471          501        278   1905.986926            NA   True

    Second example: multiple cuts at low confidence edges

    >>> # Get neuron
    >>> n = pymaid.get_neuron(27295)
    >>> # Get IDs of low confidence treenodes
    >>> lc = n.nodes[n.nodes.confidence < 5].treenode_id.values
    >>> # Cut neuron
    >>> nl = pymaid.cut_neuron(n, lc)

    See Also
    --------
    :func:`pymaid.CatmaidNeuron.prune_distal_to`
    :func:`pymaid.CatmaidNeuron.prune_proximal_to`
            ``CatmaidNeuron/List`` shorthands to this function.
    :func:`pymaid.subset_neuron`
            Returns a neuron consisting of a subset of its treenodes.

    """
    if ret not in ['proximal', 'distal', 'both']:
        raise ValueError('ret must be either "proximal", "distal" or "both"!')

    if isinstance(x, core.CatmaidNeuron):
        pass
    elif isinstance(x, core.CatmaidNeuronList):
        if x.shape[0] == 1:
            x = x[0]
        else:
            logger.error(
                '%i neurons provided. Please provide only a single neuron!' % x.shape[0])
            raise Exception(
                '%i neurons provided. Please provide only a single neuron!' % x.shape[0])
    else:
        raise TypeError('Unable to process data of type "{0}"'.format(type(x)))

    # Turn cut node into iterable
    if not utils._is_iterable(cut_node):
        cut_node = [cut_node]

    # Process cut nodes (i.e. if tag)
    cn_ids = []
    for cn in cut_node:
        # If cut_node is a tag (rather than an ID), try finding that node
        if isinstance(cn, str):
            if cn not in x.tags:
                raise ValueError(
                    '#{}: Found no treenode with tag {} - please double check!'.format(x.skeleton_id, cn))
            cn_ids += x.tags[cn]
        elif cn not in x.nodes.treenode_id.values:
            raise ValueError('No treenode with ID "{}" found.'.format(cn))
        else:
            cn_ids.append(cn)

    # Remove duplicates while retaining order - set() would mess that up
    seen = set()
    cn_ids = [cn for cn in cn_ids if not (cn in seen or seen.add(cn))]

    # Warn if not all returned
    if len(cn_ids) > 1 and ret != 'both':
        logger.warning('Multiple cuts should use `ret = True`.')

    # Go over all cut_nodes -> order matters!
    res = [x]
    for cn in cn_ids:
        # First, find out in which neuron the cut node is
        to_cut = [n for n in res if cn in n.nodes.treenode_id.values][0]
        to_cut_ix = res.index(to_cut)

        # Remove this neuron from results (will be cut into two)
        res.remove(to_cut)

        # Cut neuron
        if x.igraph and config.use_igraph:
            cut = _cut_igraph(to_cut, cn, ret)
        else:
            cut = _cut_networkx(to_cut, cn, ret)

        # If ret != 'both', we will get only a single neuron
        if not utils._is_iterable(cut):
            cut = [cut]

        # Add results back to results at same index, proximal first
        for c in cut[::-1]:
            res.insert(to_cut_ix, c)

    return core.CatmaidNeuronList(res)


def _cut_igraph(x, cut_node, ret):
    """Uses iGraph to cut a neuron."""
    # Make a copy
    g = x.igraph.copy()

    # Get vertex index
    cut_ix = g.vs.find(node_id=cut_node).index

    # Get edge to parent
    e = g.es.find(_source=cut_ix)

    # Remove edge
    g.delete_edges(e)

    # Make graph undirected -> otherwise .decompose() throws an error
    # This issue is fixed in the up-to-date branch of igraph-python
    # (which is not on PyPI O_o )
    g.to_undirected(combine_edges='first')

    # Get subgraph -> fastest way to get sets of nodes for subsetting
    a, b = g.decompose(mode='WEAK')
    # IMPORTANT: a,b are now UNDIRECTED graphs -> we must not keep using them!

    if x.root[0] in a.vs['node_id']:
        dist_graph, prox_graph = b, a
    else:
        dist_graph, prox_graph = a, b

    if ret == 'distal' or ret == 'both':
        dist = subset_neuron(x, dist_graph.vs['node_id'], clear_temp=False)

        # Change new root for dist
        # dist.nodes.loc[dist.nodes.treenode_id == cut_node, 'parent_id'] = None
        dist.nodes.loc[dist.nodes.treenode_id == cut_node, 'type'] = 'root'

        # Clear other temporary attributes
        dist._clear_temp_attr(exclude=['igraph', 'type', 'classify_nodes'])

    if ret == 'proximal' or ret == 'both':
        prox = subset_neuron(x, prox_graph.vs['node_id'] + [cut_node], clear_temp=False)

        # Change new root for dist
        prox.nodes.loc[prox.nodes.treenode_id == cut_node, 'type'] = 'end'

        # Clear other temporary attributes
        prox._clear_temp_attr(exclude=['igraph', 'type', 'classify_nodes'])

    if ret == 'both':
        return dist, prox
    elif ret == 'distal':
        return dist
    elif ret == 'proximal':
        return prox


def _cut_networkx(x, cut_node, ret):
    """Uses networkX graph to cut a neuron."""

    # Get subgraphs consisting of nodes distal to cut node
    dist_graph = nx.bfs_tree(x.graph, cut_node, reverse=True)

    if ret == 'distal' or ret == 'both':
        # bfs_tree does not preserve 'weight'
        # -> need to subset original graph by those nodes
        dist_graph = x.graph.subgraph(dist_graph.nodes)

        # Generate new neurons
        # This is the actual bottleneck of the function: ~70% of time
        dist = subset_neuron(x, dist_graph, clear_temp=False)

        # Change new root for dist
        dist.nodes.loc[dist.nodes.treenode_id == cut_node, 'parent_id'] = None
        dist.nodes.loc[dist.nodes.treenode_id == cut_node, 'type'] = 'root'

        # Reassign graphs
        dist.graph = dist_graph

        # Clear other temporary attributes
        dist._clear_temp_attr(exclude=['graph', 'type', 'classify_nodes'])

    if ret == 'proximal' or ret == 'both':
        # bfs_tree does not preserve 'weight'
        # need to subset original graph by those nodes
        prox_graph = x.graph.subgraph(
            [n for n in x.graph.nodes if n not in dist_graph.nodes] + [cut_node])

        # Generate new neurons
        # This is the actual bottleneck of the function: ~70% of time
        prox = subset_neuron(x, prox_graph, clear_temp=False)

        # Change cut node to end node for prox
        prox.nodes.loc[prox.nodes.treenode_id == cut_node, 'type'] = 'end'

        # Reassign graphs
        prox.graph = prox_graph

        # Clear other temporary attributes
        prox._clear_temp_attr(exclude=['graph', 'type', 'classify_nodes'])

    # ATTENTION: prox/dist_graph contain pointers to the original graph
    # -> changes to attributes will propagate back

    if ret == 'both':
        return dist, prox
    elif ret == 'distal':
        return dist
    elif ret == 'proximal':
        return prox


def subset_neuron(x, subset, clear_temp=True, remove_disconnected=True,
                  inplace=False):
    """ Subsets a neuron to a set of treenodes.

    Parameters
    ----------
    x :                   CatmaidNeuron
    subset :              np.ndarray | NetworkX.Graph
                          Treenodes to subset the neuron to.
    clear_temp :          bool, optional
                          If True, will reset temporary attributes (graph,
                          node classification, etc. ). In general, you should
                          leave this at ``True``.
    remove_disconnected : bool, optional
                          If True, will remove disconnected connectors that
                          have "lost" their parent treenode.
    inplace :             bool, optional
                          If False, a copy of the neuron is returned.

    Returns
    -------
    CatmaidNeuron

    Examples
    --------
    Subset neuron to presynapse-bearing branches

    >>> # Get neuron
    >>> n = pymaid.get_neuron(16)
    >>> # Go over each segment and find those with synapses
    >>> cn_nodes = set(n.presynapses.treenode_id.values)
    >>> syn_segs = [s for s in n.small_segments if set(s) & set(cn_nodes)]
    >>> # Flatten segments into list of nodes
    >>> syn_branches = [n for s in syn_segs for n in s]
    >>> # Subset neuron
    >>> axon = pymaid.subset_neuron(n, syn_branches)

    See Also
    --------
    :func:`~pymaid.cut_neuron`
            Cut neuron at specific points.

    """
    if isinstance(x, core.CatmaidNeuronList) and len(x) == 1:
        x = x[0]

    if not isinstance(x, core.CatmaidNeuron):
        raise TypeError('Can only process data of type "CatmaidNeuron", not\
                         "{0}"'.format(type(x)))

    if isinstance(subset, np.ndarray):
        pass
    elif isinstance(subset, (list, set)):
        subset = np.array(subset)
    elif isinstance(subset, (nx.DiGraph, nx.Graph)):
        subset = subset.nodes
    else:
        raise TypeError('Can only subset to list, set, numpy.ndarray or \
                         networkx.Graph, not "{0}"'.format(type(subset)))

    # Make a copy of the neuron
    if not inplace:
        x = x.copy(deepcopy=False)

    # Filter treenodes
    x.nodes = x.nodes[x.nodes.treenode_id.isin(subset)]

    # Make sure that there are root nodes
    # This is the fastest "pandorable" way: instead of overwriting the column,
    # concatenate a new column to this DataFrame
    x.nodes = pd.concat([x.nodes.drop('parent_id', axis=1),
                         x.nodes[['parent_id']].where(x.nodes.parent_id.isin(x.nodes.treenode_id.values),
                                                      None, inplace=False)],
                        axis=1)

    # Filter connectors
    if remove_disconnected:
        x.connectors = x.connectors[x.connectors.treenode_id.isin(subset)]

    # Filter tags
    x.tags = {t: [tn for tn in x.tags[t] if tn in subset] for t in x.tags}

    # Remove empty tags
    x.tags = {t: x.tags[t] for t in x.tags if x.tags[t]}

    # Fix graph representations
    if 'graph' in x.__dict__:
        x.graph = x.graph.subgraph(x.nodes.treenode_id.values)
    if 'igraph' in x.__dict__:
        if x.igraph and config.use_igraph:
            id2ix = {n: ix for ix, n in zip(x.igraph.vs.indices,
                                            x.igraph.vs.get_attribute_values('node_id'))}
            indices = [id2ix[n] for n in x.nodes.treenode_id.values]
            vs = x.igraph.vs[indices]
            x.igraph = x.igraph.subgraph(vs)

    # Reset indices of data tables
    x.nodes.reset_index(inplace=True, drop=True)
    x.connectors.reset_index(inplace=True, drop=True)

    # Clear temporary attributes
    if clear_temp:
        x._clear_temp_attr(exclude=['graph', 'igraph'])

    return x


def generate_list_of_childs(x):
    """ Returns list of childs.

    Parameters
    ----------
    x :     CatmaidNeuron | CatmaidNeuronList
            If List, must contain a SINGLE neuron.

    Returns
    -------
    dict
        ``{ treenode_id : [ child_treenode, child_treenode, ... ] }``

    """

    return {n: [e[0] for e in x.graph.in_edges(n)] for n in x.graph.nodes}


def node_label_sorting(x):
    """ Returns treenodes ordered by node label sorting according to Cuntz
    et al., PLoS Computational Biology (2010).

    Parameters
    ----------
    x :         CatmaidNeuron

    Returns
    -------
    list
        ``[ root, treenode_id, treenode_id, ... ]``

    """
    if not isinstance(x, core.CatmaidNeuron):
        raise TypeError('Need CatmaidNeuron, got "{0}"'.format(type(x)))

    if len(x.root) > 1:
        raise ValueError('Unable to process multi-root neurons!')

    # Get relevant branch points
    term = x.nodes[x.nodes.type == 'end'].treenode_id.values

    # Get distance from all branch_points
    geo = geodesic_matrix(x, tn_ids=term, directed=True)
    # Set distance between unreachable points to None
    # Need to reinitialise SparseMatrix to replace float('inf') with NaN
    #dist_mat[geo == float('inf')] = None
    dist_mat = pd.SparseDataFrame(np.where(geo == float('inf'),
                                           np.nan,
                                           geo),
                                  columns=geo.columns,
                                  index=geo.index)

    # Get starting points and sort by longest path to a terminal
    curr_points = sorted(list(x.simple.graph.predecessors(x.root[0])),
                         key=lambda n: dist_mat[n].max(),
                         reverse=True)

    # Walk from root along towards terminals, prioritising longer branches
    nodes_walked = []
    while curr_points:
        nodes_walked.append(curr_points.pop(0))
        if nodes_walked[-1] in term:
            pass
        else:
            new_points = sorted(list(x.simple.graph.predecessors(nodes_walked[-1])),
                                key=lambda n: dist_mat[n].max(),
                                reverse=True)
            curr_points = new_points + curr_points

    # Translate into segments
    node_list = [x.root[0]]
    segments = _break_segments(x)
    for n in nodes_walked:
        node_list += [seg for seg in segments if seg[0] == n][0][:-1]

    return node_list


def _igraph_to_sparse(graph, weight_attr=None):
    edges = graph.get_edgelist()
    if weight_attr is None:
        weights = [1] * len(edges)
    else:
        weights = graph.es[weight_attr]
    if not graph.is_directed():
        edges.extend([(v, u) for u, v in edges])
        weights.extend(weights)
    return csr_matrix((weights, zip(*edges)),
                      shape=(len(graph.vs), len(graph.vs)))
