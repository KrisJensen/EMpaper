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


""" This module contains functions to manipulate neuron morphology.
"""

import pandas as pd
import numpy as np
import scipy
import scipy.interpolate
from pymaid import core, graph_utils, config

# Set up logging
logger = config.logger

__all__ = sorted(['downsample_neuron', 'resample_neuron'])


def resample_neuron(x, resample_to, method='linear', inplace=False,
                    skip_errors=True):
    """ Resamples neuron(s) to given NM resolution. Preserves root, leafs,
    branchpoints. Tags and connectors are mapped onto the closest
    new treenode. Columns "confidence" and "creator" of the treenode table
    are discarded.

    Important
    ---------
    This generates an entirely new set of treenode IDs! Those will be unique
    within a neuron, but you may encounter duplicates across neurons.

    Also: be aware that high-resolution neurons will use A LOT of memory.

    Parameters
    ----------
    x :                 CatmaidNeuron | CatmaidNeuronList
                        Neuron(s) to resample.
    resample_to :       int
                        New resolution in NANOMETERS.
    method :            str, optional
                        See ``scipy.interpolate.interp1d`` for possible options.
                        By default, we're using linear interpolation.
    inplace :           bool, optional
                        If True, will modify original neuron. If False, a
                        resampled copy is returned.
    skip_errors :       bool, optional
                        If True, will skip errors during interpolation and
                        only print summary.


    Returns
    -------
    CatmaidNeuron/List
                        Downsampled neuron(s). Only if ``inplace=False``.

    See Also
    --------
    :func:`pymaid.downsample_neuron`
                        This function reduces the number of nodes instead of
                        resample to certain resolution. Useful if you are
                        just after some simplification e.g. for speeding up
                        your calculations or you want to preserve treenode IDs.
    """

    if isinstance(x, core.CatmaidNeuronList):
        results = [resample_neuron(x[i], resample_to,
                                  method='method', inplace=inplace,
                                  skip_errors=skip_errors)
                   for i in config.trange(x.shape[0],
                                          desc='Resampl. neurons',
                                          disable=config.pbar_hide,
                                          leave=config.pbar_leave)]
        if not inplace:
            return core.CatmaidNeuronList(results)
    elif not isinstance(x, core.CatmaidNeuron):
        logger.error('Unexpected datatype: %s' % str(type(x)))
        raise ValueError

    if not inplace:
        x = x.copy()

    # Collect some information for later
    nodes = x.nodes.set_index('treenode_id')
    locs = nodes[['x', 'y', 'z']]
    radii = nodes['radius'].to_dict()

    new_nodes = []
    max_tn_id = x.nodes.treenode_id.max() + 1

    errors = 0

    # Iterate over segments
    for i, seg in enumerate(config.tqdm(x.small_segments, desc='Proc. segments',
                                        disable=config.pbar_hide, leave=False)):
        # Get coordinates
        coords = locs.loc[seg].values.astype(float)
        # Get radii
        rad = [radii[tn] for tn in seg]

        # Vecs between subsequently measured points
        vecs = np.diff(coords.T)

        # path: cum distance along points (norm from first to ith point)
        path = np.cumsum(np.linalg.norm(vecs, axis=0))
        path = np.insert(path, 0, 0)

        # If path is too short, just keep the first and last treenode
        if path[-1] < resample_to or (method == 'cubic' and len(seg) <= 3):
            new_nodes += [[seg[0], seg[-1], None, coords[0]
                           [0], coords[0][1], coords[0][2], radii[seg[0]], 5]]
            continue

        # Coords of interpolation
        n_nodes = int(path[-1] / resample_to)
        interp_coords = np.linspace(path[0], path[-1], n_nodes)

        try:
            sampleX = scipy.interpolate.interp1d(path, coords[:, 0],
                                                 kind=method)
            sampleY = scipy.interpolate.interp1d(path, coords[:, 1],
                                                 kind=method)
            sampleZ = scipy.interpolate.interp1d(path, coords[:, 2],
                                                 kind=method)
            sampleR = scipy.interpolate.interp1d(path, rad,
                                                 kind=method)
        except ValueError as e:
            if skip_errors:
                errors += 1
                new_nodes += x.nodes.loc[x.nodes.treenode_id.isin(seg[:-1]),
                                         ['treenode_id', 'parent_id',
                                          'creator_id', 'x', 'y', 'z',
                                          'radius', 'confidence']].values.tolist()
                continue
            else:
                raise e

        # Sample each dim
        xnew = sampleX(interp_coords)
        ynew = sampleY(interp_coords)
        znew = sampleZ(interp_coords)
        rnew = sampleR(interp_coords).round(1)

        # Generate new coordinates
        new_coords = np.array([xnew, ynew, znew]).T.round()

        # Generate new ids (start and end node IDs of this segment)
        new_ids = seg[:1] + [max_tn_id +
                             i for i in range(len(new_coords) - 2)] + seg[-1:]

        # Keep track of new nodes
        new_nodes += [[tn, pn, None, co[0], co[1], co[2], -1, 5]
                      for tn, pn, co, r in zip(new_ids[:-1],
                                               new_ids[1:],
                                               new_coords,
                                               rnew)]

        # Increase max index
        max_tn_id += len(new_ids)

    if errors:
        logger.warning(
            '{} ({:.0%}) segments skipped due to errors'.format(errors, errors / i))

    # Add root node(s)
    root = x.root
    if not isinstance(root, (np.ndarray, list)):
        root = [x.root]
    root = x.nodes.loc[x.nodes.treenode_id.isin(root), [
        'treenode_id', 'parent_id', 'creator_id', 'x', 'y',
        'z', 'radius', 'confidence']]
    new_nodes += [list(r) for r in root.values]

    # Generate new nodes dataframe
    new_nodes = pd.DataFrame(data=new_nodes,
                             columns=['treenode_id', 'parent_id', 'creator_id',
                                      'x', 'y', 'z', 'radius', 'confidence'],
                             dtype=object
                             )

    # Convert columns to appropriate dtypes
    dtypes = {'treenode_id': int, 'parent_id': object, 'x': int, 'y': int,
              'z': int, 'radius': int, 'confidence': int}

    for k, v in dtypes.items():
        new_nodes[k] = new_nodes[k].astype(v)

    # Remove duplicate treenodes (branch points)
    new_nodes = new_nodes[~new_nodes.treenode_id.duplicated()]

    # Map connectors back:
    # 1. Get position of old synapse-bearing treenodes
    old_tn_position = x.nodes.set_index(
        'treenode_id').loc[x.connectors.treenode_id, ['x', 'y', 'z']].values
    # 2. Get closest neighbours
    distances = scipy.spatial.distance.cdist(
        old_tn_position, new_nodes[['x', 'y', 'z']].values)
    min_ix = np.argmin(distances, axis=1)
    # 3. Map back onto neuron
    x.connectors['treenode_id'] = new_nodes.iloc[min_ix].treenode_id.values

    # Map tags back:
    if x.tags:
        # 1. Get position of old tag bearing treenodes
        tag_tn = set([tn for l in x.tags.values() for tn in l])
        old_tn_position = x.nodes.set_index(
            'treenode_id').loc[tag_tn, ['x', 'y', 'z']].values
        # 2. Get closest neighbours
        distances = scipy.spatial.distance.cdist(
            old_tn_position, new_nodes[['x', 'y', 'z']].values)
        min_ix = np.argmin(distances, axis=1)
        # 3. Create a dictionary
        new_tag_tn = {
            tn: new_nodes.iloc[min_ix[i]].treenode_id for i, tn in enumerate(tag_tn)}
        # 4. Map tags back
        new_tags = {t: [new_tag_tn[tn] for tn in x.tags[t]] for t in x.tags}
        x.tags = new_tags

    # Set nodes
    x.nodes = new_nodes

    # Clear and regenerate temporary attributes
    x._clear_temp_attr()

    if not inplace:
        return x


def downsample_neuron(x, resampling_factor, preserve_cn_treenodes=True,
                      preserve_tag_treenodes=False, inplace=False,):
    """ Downsamples neuron(s) by a given factor. Preserves root, leafs,
    branchpoints by default. Preservation of treenodes with synapses can
    be toggled.

    Parameters
    ----------
    x :                      CatmaidNeuron | CatmaidNeuronList
                             Neuron(s) to downsample.
    resampling_factor :      int
                             Factor by which to reduce the node count.
    preserve_cn_treenodes :  bool, optional
                             If True, treenodes that have connectors are
                             preserved.
    preserve_tag_treenodes : bool, optional
                             If True, treenodes with tags are preserved.
    inplace :                bool, optional
                             If True, will modify original neuron. If False, a
                             downsampled copy is returned.

    Returns
    -------
    CatmaidNeuron/List
                             Downsampled neuron. Only if ``inplace=False``.

    Notes
    -----
    Use ``resampling_factor=float('inf')`` and ``preserve_cn_treenodes=False``
    to get a neuron consisting only of root, branch and end points.

    See Also
    --------
    :func:`pymaid.resample_neuron`
                             This function resamples a neuron to given
                             resolution. This will not preserve treenode IDs!

    """
    if isinstance(x, core.CatmaidNeuronList):
        return core.CatmaidNeuronList([downsample_neuron(n, resampling_factor, inplace=inplace) for n in x])
    elif isinstance(x, core.CatmaidNeuron):
        if not inplace:
            x = x.copy()
    else:
        logger.error('Unexpected datatype: %s' % str(type(x)))
        raise ValueError

    # If no resampling, simply return neuron
    if resampling_factor <= 1:
        raise ValueError('Resampling factor must be > 1.')

    if x.nodes.shape[0] <= 1:
        logger.warning(
            'No nodes in neuron {}. Skipping...'.format(x.skeleton_id))
        if not inplace:
            return x
        else:
            return

    logger.debug('Preparing to downsample neuron...')

    list_of_parents = {
        n.treenode_id: n.parent_id for n in x.nodes.itertuples()}
    list_of_parents[None] = None

    if 'type' not in x.nodes:
        graph_utils.classify_nodes(x)

    selection = x.nodes.type != 'slab'

    if preserve_cn_treenodes:
        selection = selection | x.nodes.treenode_id.isin(
            x.connectors.treenode_id)

    if preserve_tag_treenodes:
        with_tags = [t for l in x.tags.values() for t in l]
        selection = selection | x.nodes.treenode_id.isin(with_tags)

    fix_points = x.nodes[selection].treenode_id.values

    # Add soma node
    if not isinstance(x.soma, type(None)) and x.soma not in fix_points:
        fix_points = np.append(fix_points, x.soma)

    # Walk from all fix points to the root - jump N nodes on the way
    new_parents = {}

    logger.debug(
        'Sampling neuron down by factor of {0}'.format(resampling_factor))

    for en in fix_points:
        this_node = en

        while True:
            stop = False
            new_p = list_of_parents[this_node]
            if new_p:
                i = 0
                while i < resampling_factor:
                    if new_p in fix_points or not new_p:
                        new_parents[this_node] = new_p
                        stop = True
                        break
                    new_p = list_of_parents[new_p]
                    i += 1

                if stop is True:
                    break
                else:
                    new_parents[this_node] = new_p
                    this_node = new_p
            else:
                new_parents[this_node] = None
                break

    new_nodes = x.nodes[x.nodes.treenode_id.isin(
        list(new_parents.keys()))].copy()
    new_nodes.loc[:, 'parent_id'] = [new_parents[tn]
                                     for tn in new_nodes.treenode_id]

    # We have to temporarily set parent of root node from 1 to an integer
    root_ix = new_nodes[new_nodes.parent_id.isnull()].index
    new_nodes.loc[root_ix, 'parent_id'] = 0
    # first convert everything to int
    new_nodes.loc[:, 'parent_id'] = new_nodes.parent_id.values.astype(int)
    # then back to object so that we can add a 'None'
    new_nodes.loc[:, 'parent_id'] = new_nodes.parent_id.values.astype(object)

    # Reassign parent_id None to root node
    new_nodes.loc[root_ix, 'parent_id'] = None

    logger.debug(
        'Nodes before/after: {}/{}'.format(len(x.nodes), len(new_nodes)))

    x.nodes = new_nodes

    # This is essential -> otherwise e.g. graph.neuron2graph will fail
    x.nodes.reset_index(inplace=True, drop=True)

    x._clear_temp_attr()

    if not inplace:
        return x
