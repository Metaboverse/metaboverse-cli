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
import networkx as nx
import pandas as pd
import requests
import pickle
import json
import sys
import os

"""Import internal dependencies
"""
try:
    from analyze.prepare_data import __main__ as prepare_data
    from analyze.model import __template__
    from analyze.model import __model__
    from analyze.model import load_references
    from analyze.model import load_metabolite_synonym_dictionary
    from utils import progress_feed, read_network, get_metaboverse_cli_version
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
    prepare_data = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prepare_data)
    prepare_data = prepare_data.__main__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/analyze/model.py"))
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    __template__ = model.__template__
    __model__ = model.__model__
    load_references = model.load_references
    load_metabolite_synonym_dictionary = model.load_metabolite_synonym_dictionary

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    read_network = utils.read_network
    get_metaboverse_cli_version = utils.get_metaboverse_cli_version


SOURCEFORGE_URL='https://sourceforge.net/projects/metaboverse/files/mvrs_files/'


def process_data(
        network,
        args_dict):
    """
    """
    if str(args_dict['transcriptomics']).lower() != 'none' \
            or str(args_dict['proteomics']).lower() != 'none' \
            or str(args_dict['metabolomics']).lower() != 'none':

        data, stats, unmapped = prepare_data(
            network=network,
            transcriptomics_url=args_dict['transcriptomics'],
            proteomics_url=args_dict['proteomics'],
            metabolomics_url=args_dict['metabolomics'],
            database_source=args_dict['database_source'])
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
        flag_data = True

    return data, stats, unmapped, flag_data


def read_template(
        args_dict,
        network,
        url):
    """
    """
    print('Downloading Metaboverse graph template for organism...')

    file = os.path.join(
        args_dict['output'],
        args_dict['organism_id'] + '_template.mvrs')
    os.system('curl -L ' + url + ' -o \"' + file + '\"')

    with open(file) as graph_template:
        graph_data = json.load(graph_template)

    graph = nx.readwrite.json_graph.node_link_graph(
        {
            'nodes': graph_data['nodes'],
            'links': graph_data['links']
        },
        directed=graph_data['directed'],
        multigraph=graph_data['multigraph'])
    network['reaction_database'] = graph_data['reaction_dictionary']
    network['pathway_database'] = graph_data['pathway_dictionary']
    super_pathways = graph_data['super_pathways']
    degree_dictionary = graph_data['degree_dictionary']

    reverse_genes, protein_dictionary, chebi_dictionary, \
    name_reference, uniprot_mapper = load_references(
        args_dict=args_dict,
        ensembl=network['ensembl_synonyms'],
        uniprot=network['uniprot_synonyms'],
        chebi=network['chebi_mapper'],
        uniprot_metabolites=network['uniprot_metabolites'])
    metabolite_mapper = load_metabolite_synonym_dictionary()

    progress_feed(args_dict, "model", 9)

    return graph, args_dict, network, name_reference, \
        degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper


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
    data, stats, unmapped, flag_data = process_data(
        network=network,
        args_dict=args_dict)
    progress_feed(args_dict, "model", 3)

    # Generate graph template
    this_version = get_metaboverse_cli_version()
    test_url = (
        SOURCEFORGE_URL
        + this_version + '/'
        + args_dict['organism_id'] + '_template.mvrs/download')
    url_response = requests.head(test_url)

    if (args_dict['force_new_curation'] == False \
    or args_dict['force_new_curation'] == "False") \
    and url_response.status_code != 404:
        graph, args_dict, network, name_reference, \
        degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper = read_template(
            args_dict=args_dict,
            network=network,
            url=test_url)
    else:
        graph, args_dict, network, name_reference, \
        degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper = __template__(
            args_dict=args_dict,
            network=network,
            species_id=args_dict['organism_id'],
            output_file=args_dict['output_file'])

    # Overlay data on graph and collapse as able
    graph_name = __model__(
        graph=graph,
        args_dict=args_dict,
        network=network,
        data=data,
        stats=stats,
        species_id=args_dict['organism_id'],
        output_file=args_dict['output_file'],
        name_reference=name_reference,
        degree_dictionary=degree_dictionary,
        chebi_dictionary=chebi_dictionary,
        uniprot_mapper=uniprot_mapper,
        metabolite_mapper=metabolite_mapper,
        super_pathways=super_pathways,
        unmapped=unmapped,
        flag_data=flag_data)

    progress_feed(args_dict, "model", 10)


def test():
    network = read_network(
        network_url="C:\\Users\\jorda\\Desktop\\HSA.mvdb")
    len(list(network['reaction_database'].keys()))

    args_dict = {
        'output': "C:\\Users\\jorda\\Desktop",
        'organism_id': 'HSA'}


    graph
    adj_matrix = nx.linalg.graphmatrix.adjacency_matrix(graph).todense()

    len(list(graph.nodes()))

    type(adj_matrix[0])

    import pandas as pd
    df = pd.DataFrame(
            adj_matrix,
            index=list(graph.nodes()),
            columns=list(graph.nodes())).apply(pd.to_numeric)

    # make this dictionary with the template stage to save a couple minutes of processing time
    df_copy = df.copy()
    #df_copy.iloc[:, :].replace(1, pd.Series(df_copy.columns, df_copy.columns))

    neighbor_dict = {}
    col_labels = df_copy.columns.tolist()

    for name, row in df_copy.iterrows():
        indices = [i for i, x in enumerate(row) if x == 1]
        neighbor_dict[name] = [col_labels[_i] for _i in indices]
