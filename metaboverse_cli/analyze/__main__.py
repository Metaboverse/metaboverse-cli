"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

Copyright (C) 2019-2021 Jordan A. Berg
Email: jordan<dot>berg<at>biochem<dot>utah<dot>edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function
import pandas as pd
import pickle

"""Import internal dependencies
"""
try:
    from analyze.prepare_data import __main__ as prepare_data
    from analyze.model import __main__ as model
    from utils import progress_feed, read_network
except:
    import os
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
    prepare_data = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prepare_data)
    prepare_data = prepare_data.__main__

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/analyze/model.py"))
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    model = model.__main__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    read_network = utils.read_network


def __main__(
        args_dict):
    """Analyze data on network model
    """

    # Get network curation info
    network = read_network(
        network_url=args_dict['network'])
    progress_feed(args_dict, "model", 2)

    if args_dict['organism_curation'] != 'None':
        args_dict['organism_id'] = network['organism_id']

    # Read in data (if any)
    if str(args_dict['transcriptomics']).lower() != 'none' \
            or str(args_dict['proteomics']).lower() != 'none' \
            or str(args_dict['metabolomics']).lower() != 'none':

        data, stats, unmapped = prepare_data(
            network=network,
            transcriptomics_url=args_dict['transcriptomics'],
            proteomics_url=args_dict['proteomics'],
            metabolomics_url=args_dict['metabolomics'],
            database_source=args_dict['database_source'])
        progress_feed(args_dict, "model", 3)
        flag_data = False

    else:
        data = pd.DataFrame()
        data['NoSample'] = [0, 0, 0]
        data.index = ['dummy_index1', 'dummy_index2', 'dummy_index3']

        stats = pd.DataFrame()
        stats['NoSample'] = [0, 0, 0]
        stats.index = ['dummy_index1', 'dummy_index2', 'dummy_index3']

        unmapped = {
            'transcriptomics_unmapped': [],
            'proteomics_unmapped': []
        }

        progress_feed(args_dict, "model", 3)
        flag_data = True

    # Generate graph
    graph_name = model(
        args_dict=args_dict,
        network=network,
        data=data,
        stats=stats,
        species_id=args_dict['organism_id'],
        output_file=args_dict['output_file'],
        unmapped=unmapped,
        flag_data=flag_data)
    progress_feed(args_dict, "model", 10)
