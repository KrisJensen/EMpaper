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

""" This module contains functions to save and restore neuron data.
"""

from pymaid import core, fetch, utils, config, connectivity
import pandas as pd
import datetime

# Set up logging
logger = config.logger

__all__ = ['take_snapshot', 'load_snapshot']


def take_snapshot(x, skeleton_data=True, cn_table=False, node_details=False,
                  adjacency_matrix=True, remote_instance=None,
                  cn_details=True, annotations=False):
    """ Take a snapshot of CATMAID data associated with a set of neurons.

    Important
    ---------
    If you pass Catmaidneuron/List that have been modified (e.g. pruned),
    other data (e.g. connectivity, etc) will be subset as well if applicable.
    If your CatmaidNeuron/List is still naive, you might want to just pass the
    skeleton ID(s) to speed things up.

    Parameters
    ----------
    x :                 skeleton IDs | CatmaidNeuron/List
                        Neurons for which to retrieve data.
    skeleton_data :     bool, optional
                        Include 3D skeleton data.
    cn_table :          bool, optional
                        Include connectivity table. Covers all neurons
                        connected to input neurons.
    node_details :      bool, optional
                        Include treenode and connector details.
    adjacency_matrix :  bool, optional
                        Include adjacency matrix covering the input neurons.
    cn_details :        bool, optional
                        Include connector details.
    annotations :       bool, optional
                        Include neuron annotations.
    remote_instance :   Catmaid Instance, optional
                        Either pass explicitly or define globally. Will
                        obviously not be added to the snapshot!

    Returns
    -------
    pandas Series

    Examples
    --------

    See Also
    --------
    :func:`~pymaid.load_snapshot`
            Use to load a snapshot file.
    """

    remote_instance = utils._eval_remote_instance(remote_instance)

    skids = utils.eval_skids(x, remote_instance)

    snapshot = pd.Series()

    # Add Coordinates Universal Time
    snapshot['utc_date'] = datetime.datetime.utcnow()

    # Add pymaid version
    snapshot['pymaid_version'] = init.__version__

    # Add skeleton data
    if skeleton_data:
        if not isinstance(x, (core.CatmaidNeuronList, core.CatmaidNeuron)):
            skdata = fetch.get_neurons(skids, remote_instance=remote_instance)
        else:
            skdata = x

        if isinstance(skdata, core.CatmaidNeuron):
            skdata = core.CatmaidNeuronList(skdata)

        snapshot['skeleton_data'] = skdata.to_dataframe()

    # Add connectivity table
    if cn_table:
        if isinstance(x, (core.CatmaidNeuron, core.CatmaidNeuronList)):
            snapshot['cn_table'] = connectivity.cn_table_from_connectors(x, remote_instance=remote_instance)
        else:
            # Add connectivity table
            snapshot['cn_table'] = fetch.get_partners(x, remote_instance=remote_instance)

    # Add connectivity table
    if node_details:
        snapshot['treenode_details'] = fetch.get_node_details(skdata.nodes.treenode_id.values,
                                                              remote_instance=remote_instance)
        snapshot['connector_details'] = fetch.get_node_details(skdata.connectors.connector_id.values,
                                                               remote_instance=remote_instance)

    # Add adjacency matrix
    if adjacency_matrix:
        if isinstance(x, (core.CatmaidNeuron, core.CatmaidNeuronList)):
            snapshot['adjacency_matrix'] = connectivity.adjacency_from_connectors(x, remote_instance=remote_instance)
        else:
            # Add connectivity table
            snapshot['adjacency_matrix'] = connectivity.adjacency_matrix(x, remote_instance=remote_instance)

    # Add annotations
    if annotations:
        snapshot['annotations'] = fetch.get_annotations(skids)

    # Add connector details
    if cn_details:
        snapshot['cn_details'] = fetch.get_connector_details(skdata)

    return snapshot


def load_snapshot(x):
    pass
