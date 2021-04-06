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
    from analyze.utils import remove_defective_reactions
    from utils import progress_feed, track_progress, read_network, \
                      get_metaboverse_cli_version, write_database, safestr
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
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

    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "analyze", "utils.py"
                     ))
    spec = importlib.util.spec_from_file_location("", module_path)
    analyze_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyze_utils)
    remove_defective_reactions = analyze_utils.remove_defective_reactions

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    track_progress = utils.track_progress
    read_network = utils.read_network
    get_metaboverse_cli_version = utils.get_metaboverse_cli_version
    write_database = utils.write_database
    safestr = utils.safestr


TEMPLATE_URL='https://sourceforge.net/projects/metaboverse/files/mvrs_files/'
NEIGHBOR_URL='https://sourceforge.net/projects/metaboverse/files/nbdb_files/'


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
        url,
        user_provided=False):
    """
    """
    def get_template_file(url, args_dict):
        file = os.path.join(
            args_dict['output'],
            args_dict['organism_id'] + '_template.mvrs')
        os.system('curl -L ' + url + ' -o \"' + file + '\"')
        return file

    print('Downloading Metaboverse graph template for organism...')

    if user_provided == False:
        file = get_template_file(url, args_dict)
    else:
        file = url

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

    progress_feed(args_dict, "graph", 9)

    return graph, args_dict, network, name_reference, \
        degree_dictionary, super_pathways, chebi_dictionary, \
        uniprot_mapper, metabolite_mapper


def download_neighbors_dictionary(
        args_dict,
        url,
        user_provided=False):
    """
    """
    print('Downloading Metaboverse neighbors dictionary for organism...')

    def get_neighbor_file(url, args_dict):
        file = os.path.join(
            args_dict['output'],
            args_dict['organism_id'] + '.nbdb')
        os.system('curl -L ' + url + ' -o \"' + file + '\"')
        return file

    if user_provided == False:
        file = get_neighbor_file(url, args_dict)
    else:
        file = url

    neighbors_dictionary = read_network(
        file_path=args_dict['output'],
        network_url=file)

    return neighbors_dictionary


def make_neighbors_dictionary(
        args_dict,
        graph,
        reaction_dictionary):
    """
    """
    reaction_ids = set(reaction_dictionary.keys())

    print('Generating Metaboverse neighbors dictionary for organism...')

    counter = 0
    edges = list(graph.edges)
    edge_len = len(edges)
    neighbors_dictionary = {}
    for e in edges:
        one = e[0]
        two = e[1]
        if one in neighbors_dictionary.keys():
            neighbors_dictionary[one].add(two)
        else:
            neighbors_dictionary[one] = set()
            neighbors_dictionary[one].add(two)
        if two in neighbors_dictionary.keys():
            neighbors_dictionary[two].add(one)
        else:
            neighbors_dictionary[two] = set()
            neighbors_dictionary[two].add(one)
        counter = track_progress(args_dict, counter, edge_len, 5)

    print('Tuning neighbors dictionary...')
    reaction_neighbors_dictionary = {}
    counter = 0
    neighbors = list(neighbors_dictionary.keys())
    neighbors_number = len(neighbors)
    for neighbor in neighbors:
        if neighbor in reaction_ids:
            components = neighbors_dictionary[neighbor]
            connected_reactions = set()
            for _c in components:
                for _c_ in neighbors_dictionary[_c]:
                    if _c_ in reaction_ids:
                        connected_reactions.add(_c_)
            connected_reactions = list(connected_reactions)
            reaction_neighbors_dictionary[neighbor] = connected_reactions
        counter = track_progress(args_dict, counter, neighbors_number, 5)

    print('Writing neighbors dictionary to database file...')
    write_database(
        output=args_dict['output'],
        file=args_dict['organism_id'] + '.nbdb',
        database=reaction_neighbors_dictionary)

    return reaction_neighbors_dictionary


def __main__(
        args_dict):
    """Analyze data on network model
    """

    # Get network curation info
    network = read_network(
        file_path=args_dict['output'],
        network_url=args_dict['curation'])
    progress_feed(args_dict, "graph", 1)

    if args_dict['organism_curation_file'] != 'None':
        args_dict['organism_id'] = network['organism_id']

    # Read in data (if any)
    data, stats, unmapped, flag_data = process_data(
        network=network,
        args_dict=args_dict)
    progress_feed(args_dict, "graph", 3)

    # Generate graph template
    this_version = get_metaboverse_cli_version()
    test_url = (
        TEMPLATE_URL
        + this_version + '/'
        + args_dict['organism_id'] + '_template.mvrs/download')
    url_response = requests.head(test_url)

    if (args_dict['force_new_curation'] == False \
    or args_dict['force_new_curation'] == "False") \
    and 'graph_template_file' in args_dict \
    and safestr(args_dict['graph_template_file']) != None \
    and safestr(args_dict['graph_template_file']) != 'None':
        try:
            graph, args_dict, network, name_reference, \
            degree_dictionary, super_pathways, chebi_dictionary, \
            uniprot_mapper, metabolite_mapper = read_template(
                args_dict=args_dict,
                network=network,
                url=args_dict['graph_template_file'],
                user_provided=True)
        except:
            graph, args_dict, network, name_reference, \
            degree_dictionary, super_pathways, chebi_dictionary, \
            uniprot_mapper, metabolite_mapper = __template__(
                args_dict=args_dict,
                network=network,
                species_id=args_dict['organism_id'],
                output_file=args_dict['output_file'])
    elif (args_dict['force_new_curation'] == False \
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

    if len(graph.nodes) == 0 or len(graph.edges) == 0:
        raise Exception("Unable to generate a reaction-based network based on the input organism template.")

    # Generate graph template
    neighbors_url = (
        NEIGHBOR_URL
        + this_version + '/'
        + args_dict['organism_id'] + '.nbdb/download')
    neighbor_response = requests.head(neighbors_url)

    force_neighbors = False
    if (args_dict['force_new_curation'] == False \
    or args_dict['force_new_curation'] == "False") \
    and 'neighbor_dictionary_file' in args_dict \
    and safestr(args_dict['neighbor_dictionary_file']) != None \
    and safestr(args_dict['neighbor_dictionary_file']) != 'None':
        try:
            neighbors_dictionary = download_neighbors_dictionary(
                args_dict=args_dict,
                url=args_dict['neighbor_dictionary_file'],
                user_provided=True)
        except:
            force_neighbors = True

    elif (args_dict['force_new_curation'] == False \
    or args_dict['force_new_curation'] == "False") \
    and neighbor_response.status_code != 404:
        try:
            neighbors_dictionary = download_neighbors_dictionary(
                args_dict=args_dict,
                url=neighbors_url)
        except:
            force_neighbors = True
    else:
        force_neighbors = True

    if force_neighbors == True:
        no_defective_reactions = remove_defective_reactions(
            network=network)
        neighbors_dictionary = make_neighbors_dictionary(
            args_dict=args_dict,
            graph=graph,
            reaction_dictionary=no_defective_reactions)

    # Overlay data on graph and collapse as able
    graph_name = __model__(
        graph=graph,
        args_dict=args_dict,
        network=network,
        data=data,
        stats=stats,
        species_id=args_dict['organism_id'],
        output_file=args_dict['output_file'],
        neighbors_dictionary=neighbors_dictionary,
        name_reference=name_reference,
        degree_dictionary=degree_dictionary,
        chebi_dictionary=chebi_dictionary,
        uniprot_mapper=uniprot_mapper,
        metabolite_mapper=metabolite_mapper,
        super_pathways=super_pathways,
        unmapped=unmapped,
        flag_data=flag_data)

    progress_feed(args_dict, "graph", 10)
    return graph_name


def test():
    args_dict = {
        'output': "C:\\Users\\jorda\\Desktop",
        'curation': "MODEL1604210000.mvdb",
        'organism_id': 'MODEL1604210000',
        'output_file': "C:\\Users\\jorda\\Desktop\\MODEL1604210000.mvrs",
        'labels': "0",
        'blocklist': "H+",
        'database_date': "",
        'curation_date': ""}
    args_dict = {
        'output': "C:\\Users\\jorda\\Desktop",
        'curation': "HSA.mvdb",
        'graph_template_file': "C:\\Users\\jorda\\Desktop\\HSA_template.mvrs",
        'organism_id': 'HSA'}

    network = read_network(
        file_path="C:\\Users\\jorda\\Desktop",
        network_url="MODEL1604210000.mvdb")

    neighbors = read_network(
        file_path="C:\\Users\\jorda\\Desktop",
        network_url="MODEL1604210000.nbdb")

    neighbors_dictionary = download_neighbors_dictionary(
        args_dict=args_dict,
        url=args_dict['neighbor_dictionary_file'],
        user_provided=True)
