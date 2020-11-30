"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

Copyright (C) 2019-2020 Jordan A. Berg
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

"""Import dependencies
"""
import pickle
import pandas as pd

"""Import internal dependencies
"""
try:
    from analyze.prepare_data import __main__ as prepare_data
    from analyze.model import __main__ as model
    from utils import progress_feed
except:
    import os
    import importlib.util
    spec = importlib.util.spec_from_file_location("__main__", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
    prepare_data = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prepare_data)
    prepare_data = prepare_data.__main__

    spec = importlib.util.spec_from_file_location("__main__", os.path.abspath("./metaboverse_cli/analyze/model.py"))
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    model = model.__main__

    spec = importlib.util.spec_from_file_location("progress_feed", os.path.abspath("./metaboverse_cli/utils.py"))
    progress_feed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(progress_feed)
    progress_feed = progress_feed.progress_feed

def test():

    args_dict = {
        'network': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE_metaboverse_db.pickle',
        'metabolomics': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/metabolomics_mct1_030min.txt',
        'transcriptomics': 'None',
        'proteomics': 'None',
        'organism_curation': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE_metaboverse_db.pickle',
        'organism_id': 'SCE',
        'output_file': '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE.mvrs',
    }

    __main__(
        args_dict=args_dict)

def test2():

    args_dict = {
        'network': 'C:\\Users\\jorda\\Desktop\\BMID000000141967.mvdb',
        'metabolomics': 'None',
        'transcriptomics': 'None',
        'proteomics': 'None',
        'organism_curation': 'C:\\Users\\jorda\\Desktop\\BMID000000141967.mvdb',
        'organism_id': 'BMID000000141967',
        'output_file': 'C:\\Users\\jorda\\Desktop\\BMID000000141967.mvrs',
        'collapse_with_modifiers': False,
        'labels': '0',
        'blocklist': 'H+, H2O, CO2, e-, Mn2+, K+, Na+, Mg2+, O2, H2O2'
    }

    with open(args_dict['network'], 'rb') as network_file:
        network = pickle.load(network_file)

    __main__(
        args_dict=args_dict)

def test3():

    # Check that CHEBI values map
    args_dict = {
        'database_source': 'reactome',
        'network': 'C:\\Users\\jorda\\Desktop\\HSA_metaboverse_db.pickle',
        'metabolomics': 'C:\\Users\\jorda\\Desktop\\targetedMetabolites_byChEBIs_rawAbundances_bySystem_unmapped.txt',
        'transcriptomics': 'None',
        'proteomics': 'None',
        'organism_curation': 'C:\\Users\\jorda\\Desktop\\HSA_metaboverse_db.pickle',
        'organism_id': 'HSA',
        'output_file': 'C:\\Users\\jorda\\Desktop\\HSA_chebi_mapping',
        'collapse_with_modifiers': False,
        'labels': 'a,b,c,d',
        'blocklist': 'H+, H2O, CO2, e-, Mn2+, K+, Na+, Mg2+, O2, H2O2',
        'database_version': 'test'
    }

    with open(args_dict['network'], 'rb') as network_file:
        network = pickle.load(network_file)

    __main__(
        args_dict=args_dict)

def test4():

    args_dict = {
        'database_source': 'reactome',
        'network': 'C:\\Users\\jorda\\Desktop\\HSA_metaboverse_db.pickle',
        'metabolomics': 'C:\\Users\\jorda\\Desktop\\projects\\manuscript\\data\\lung_tumor_pr000305_st000390\\lung_tumor_vs_normal_measurements.txt',
        'transcriptomics': 'None',
        'proteomics': 'None',
        'organism_curation': 'C:\\Users\\jorda\\Desktop\\HSA_metaboverse_db.pickle',
        'organism_id': 'HSA',
        'output_file': 'C:\\Users\\jorda\\Desktop\\lung_broadcasting.mvrs',
        'collapse_with_modifiers': False,
        'labels': 'a,b,c,d',
        'blocklist': 'H+, H2O, CO2, e-, Mn2+, K+, Na+, Mg2+, O2, H2O2',
        'database_version': 'test',
        'broadcast_genes': True,
        'broadcast_metabolites': True
    }

    with open(args_dict['network'], 'rb') as network_file:
        network = pickle.load(network_file)

    __main__(
        args_dict=args_dict)

def test_check_db():

    network_url = '/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_networks/SCE.mvdb'
    with open(network_url, 'rb') as network_file:
        network = pickle.load(network_file)

def read_network(
        network_url):
    """Read in network from previous curation module
    - was provided as a URL to the file and saved to args_dict['network'] in
    "curate" sub-module
    """

    with open(network_url, 'rb') as network_file:
        network = pickle.load(network_file)

    return network

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
        data['NoSample'] = [0,0,0]
        data.index = ['dummy_index1', 'dummy_index2', 'dummy_index3']

        stats = pd.DataFrame()
        stats['NoSample'] = [0,0,0]
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
