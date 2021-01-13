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

"""Import dependencies
"""
import os
import re
import requests
import zipfile
from datetime import date
import pandas as pd
import numpy as np
import json
import pickle
import itertools
import networkx as nx
from networkx.readwrite import json_graph
from collections import Counter

try:
    from analyze.model import build_chebi_reference, build_name_reference, load_metabolite_synonym_dictionary, uniprot_ensembl_reference, gather_synonyms, name_graph, compile_node_degrees
    from utils import progress_feed
except:
    import importlib.util
    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "analyze", "model.py"
    ))
    spec = importlib.util.spec_from_file_location("", module_path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    build_chebi_reference = model.build_chebi_reference
    build_name_reference = model.build_name_reference
    load_metabolite_synonym_dictionary = model.load_metabolite_synonym_dictionary
    uniprot_ensembl_reference = model.uniprot_ensembl_reference
    gather_synonyms = model.gather_synonyms
    name_graph = model.name_graph
    compile_node_degrees = model.compile_node_degrees

    module_path = os.path.abspath(
        os.path.join(".", "metaboverse_cli", "utils.py"
    ))
    spec = importlib.util.spec_from_file_location("progress_feed", module_path)
    progress_feed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(progress_feed)
    progress_feed = progress_feed.progress_feed

def get_metabolites(
        data,
        columns):
    """Make metabolite dictionary with the available synonyms from the MIDAS
        database
    """

    common_key = columns['Common_metabolite_name']
    common_names = set()
    for index, row in data.iterrows():
        common_names.add(row[common_key])

    return list(common_names)

def fetch_chebi_mappers(
        metabolites=[
            "1,3-Diaminopropane",
            "2-Ketobutyric acid",
            "2-Hydroxybutyric acid",
            "2-Methoxyestrone"],
        metaboanalyst_api="http://api.xialab.ca/mapcompounds"
        ):
    """Adapted from R from MetaboAnalyst 5.0
    - https://www.metaboanalyst.ca/docs/RTutorial.xhtml#3.2%20Compound%20name%20mapping
    - https://www.metaboanalyst.ca/docs/APIs.xhtml
    """

    input_str = "{\n\t\"queryList\": \""
    for v in metabolites:
        input_str = input_str + v + ";"
    input_str = input_str + "\",\n\t\"inputType\": \"name\"\n}"

    query_results = requests.post(
        metaboanalyst_api,
        data=input_str,
        headers={
            'Content-Type': "application/json",
            'cache-control': "no-cache",
            }
        )

    return query_results.json()

def define_mapper(
        mapper,
        data,
        columns):
    """Make metabolite dictionary with the available synonyms from the MIDAS
        database
    """

    _mapper = {}
    for x in mapper:
        _mapper[x['query']] = x

    name_key = columns['Metabolite']
    common_key = columns['Common_metabolite_name']
    kegg_key = columns['KEGG_ID']
    hmdb_key = columns['HMDB_ID']
    midas_key = columns['MIDAS_ID']
    smiles_key = columns['SMILES_metabolite']

    metabolites = {}
    for index, row in data.iterrows():
        _name = str(row[name_key])
        _common = str(row[common_key])
        _kegg = str(row[kegg_key])
        _hmdb = str(row[hmdb_key])
        _midas = str(row[midas_key])
        _smiles = str(row[smiles_key])

        if _name in metabolites.keys():
            metabolites[_kegg]['hmdb_id'].add(_hmdb)
            metabolites[_kegg]['name'].add(_name)
            metabolites[_kegg]['common_name'].add(_common)
            metabolites[_kegg]['kegg_id'].add(_kegg)
            metabolites[_kegg]['midas_id'].add(_midas)
            metabolites[_kegg]['smiles'].add(_smiles)
        else:
            metabolites[_kegg] = {
                'hmdb_id': set([_hmdb]),
                'name': set([_name]),
                'common_name': set([_common]),
                'kegg_id': set([_kegg]),
                'midas_id': set([_midas]),
                'smiles': set([_smiles]),
                'chebi_id': set()
            }

        if _common in _mapper.keys():
            _hmdb_mod = _mapper[_common]['hmdb_id']
            _kegg_mod = _mapper[_common]['kegg_id']
            _common_mod = _mapper[_common]['hit']
            _chebi = _mapper[_common]['chebi_id']

            metabolites[_kegg]['hmdb_id'].add(_hmdb_mod)
            metabolites[_kegg]['kegg_id'].add(_kegg_mod)
            metabolites[_kegg]['common_name'].add(_common_mod)
            metabolites[_kegg]['chebi_id'].add(_chebi)
            metabolites[_kegg]['chebi_id'].add('CHEBI:' + _chebi)

    return metabolites

def finalize_node(
        graph,
        id):
    for k in graph.nodes()[id]:
        if isinstance(graph.nodes()[id][k], set):
            graph.nodes()[id][k] = list(graph.nodes()[id][k])

    return graph

def targeted_graph(
        metabolites,
        reactions,
        pathways,
        species_reference,
        name_database,
        metabolite_mapper,
        uniprot_mapper,
        component_database):
    """Build graph
    - Add nodes and edges
    - Map names to objects in the graph for display
    - Calculate degree for each node
    """

    # Initialize graph object
    graph = nx.DiGraph()
    for _m in metabolites.keys():

        # Find species_id
        species_ids, parsed_syns = get_species(
            metabolite=metabolites[_m],
            species_reference=species_reference,
            name_database=name_database,
            metabolite_mapper=metabolite_mapper,
            uniprot_mapper=uniprot_mapper)

        # Add node
        graph.add_node(_m)
        graph.nodes()[_m]['id'] = _m
        graph.nodes()[_m]['map_id'] = _m
        graph.nodes()[_m]['name'] = metabolites[_m]['name']
        graph.nodes()[_m]['common_name'] = metabolites[_m]['common_name']
        graph.nodes()[_m]['type'] = 'metabolite'
        graph.nodes()[_m]['sub_type'] = 'midas_metabolite'
        graph.nodes()[_m]['species_ids'] = species_ids
        graph.nodes()[_m]['synonyms'] = parsed_syns
        graph.nodes()[_m]['notes'] = ""
        graph.nodes()[_m]['smiles'] = metabolites[_m]['smiles']
        graph.nodes()[_m]['pathways'] = ""

        # Find nearest neighbor reactions
        _reaction_list = set()
        for _s in species_ids:
            for _k, _v in reactions.items():
                if _s in _v['reactants'] \
                or _s in _v['products']:
                    _reaction_list.add(_k)
                for _mod in _v['modifiers']:
                    if _mod[0] == _s:
                        _reaction_list.add(_k)
        graph.nodes()[_m]['reactome_reactions'] = list(_reaction_list)
        graph = finalize_node(
            graph=graph,
            id=_m)

        # Add reactions and components
        for _r in _reaction_list:

            # Add reaction node
            graph.add_node(_r)
            graph.nodes()[_r]['id'] = _r
            graph.nodes()[_r]['map_id'] = _r
            graph.nodes()[_r]['name'] = reactions[_r]['name']
            graph.nodes()[_r]['common_name'] = ""
            graph.nodes()[_r]['type'] = 'reaction'
            graph.nodes()[_r]['sub_type'] = 'reaction'
            graph.nodes()[_r]['species_ids'] = ""
            graph.nodes()[_r]['synonyms'] = ""
            graph.nodes()[_r]['notes'] = reactions[_r]['notes']
            graph.nodes()[_m]['smiles'] = ""
            graph.nodes()[_r]['pathways'] = set()
            graph.nodes()[_r]['reactome_reactions'] = set()

            for _p in pathways.keys():
                if _r in pathways[_p]['reactions']:
                    graph.nodes()[_r]['pathways'].add(_p)
            graph.nodes()[_r]['pathways'] = list(graph.nodes()[_r]['pathways'])
            graph = finalize_node(
                graph=graph,
                id=_r)

            # For all reactant, product, and mods, add node and edge to reaction node
            _reactants = reactions[_r]['reactants']
            _products = reactions[_r]['products']
            _modifiers = reactions[_r]['modifiers']

            for _reactant in _reactants:
                _name = component_database[_reactant]['name']
                _type = 'reactant'
                _subtype = component_database[_reactant]['type']
                graph = add_node(
                    graph=graph,
                    identifier=_reactant,
                    name=_name,
                    type=_type,
                    sub_type=_subtype)
                graph = add_edge(
                    graph=graph,
                    source=_reactant,
                    target=_r,
                    type=_type,
                    sub_type=_subtype)
                graph = finalize_node(
                    graph=graph,
                    id=_reactant)

            for _product in _products:
                _name = component_database[_product]['name']
                _type = 'product'
                _subtype = component_database[_product]['type']
                graph = add_node(
                    graph=graph,
                    identifier=_product,
                    name=_name,
                    type=_type,
                    sub_type=_subtype)
                graph = add_edge(
                    graph=graph,
                    source=_r,
                    target=_product,
                    type=_type,
                    sub_type=_subtype)
                graph = finalize_node(
                    graph=graph,
                    id=_product)

            for _modifier in _modifiers:
                _name = component_database[_modifier[0]]['name']
                _type = 'modifier'
                _subtype = _modifier[1]
                graph = add_node(
                    graph=graph,
                    identifier=_modifier[0],
                    name=_name,
                    type=_type,
                    sub_type=_subtype)
                graph = add_edge(
                    graph=graph,
                    source=_modifier[0],
                    target=_r,
                    type=_type,
                    sub_type=_subtype)
                graph = finalize_node(
                    graph=graph,
                    id=_modifier[0])

    return graph

def add_edge(
        graph,
        source,
        target,
        type,
        sub_type):

    graph.add_edges_from([(source, target)])
    graph.edges()[(source, target)]['type'] = type
    graph.edges()[(source, target)]['sub_type'] = sub_type

    return graph

def add_node(
        graph,
        identifier,
        name,
        type,
        sub_type):

    graph.add_node(identifier)
    graph.nodes()[identifier]['id'] = identifier
    graph.nodes()[identifier]['map_id'] = identifier
    graph.nodes()[identifier]['name'] = name
    graph.nodes()[identifier]['common_name'] = name
    graph.nodes()[identifier]['type'] = type
    graph.nodes()[identifier]['sub_type'] = sub_type
    graph.nodes()[identifier]['species_ids'] = identifier
    graph.nodes()[identifier]['synonyms'] = ""
    graph.nodes()[identifier]['notes'] = ""
    graph.nodes()[identifier]['smiles'] = ""
    graph.nodes()[identifier]['pathways'] = set()
    graph.nodes()[identifier]['reactome_reactions'] = set()

    return graph

def get_species(
        metabolite,
        species_reference,
        name_database,
        metabolite_mapper,
        uniprot_mapper):
    """
    """

    species_ids = set()

    chebi_id = metabolite['chebi_id']
    hmdb_id = metabolite['hmdb_id']
    name = metabolite['name']
    common_name = metabolite['common_name']
    kegg_id = metabolite['kegg_id']
    identifiers = itertools.chain(
        chebi_id,
        hmdb_id,
        name,
        common_name,
        kegg_id)

    mapper_id, parsed_syns_list = gather_synonyms(
        map_id=next(iter(metabolite['name'])),
        init_syns=identifiers,
        metabolite_mapper=metabolite_mapper,
        uniprot_mapper=uniprot_mapper,
        ignore_enantiomers=True)

    # cross-reference HMDB IDs in relevant databases
    for _id in parsed_syns_list:

        if _id in species_reference:
            species_ids.add(name_database[_id])
        elif _id in name_database:
            species_ids.add(name_database[_id])

    return species_ids, parsed_syns_list

def test():

    import os
    import importlib.util
    spec = importlib.util.spec_from_file_location("", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    read_network = utils.read_network
    network = read_network(
        network_url='C:\\Users\\jorda\\Desktop\\HSA.mvdb')

    spec = importlib.util.spec_from_file_location("", os.path.abspath("./metaboverse_cli/target/utils.py"))
    build_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build_utils)
    import_midas = build_utils.import_midas
    data, columns = import_midas(
        filename='C:\\Users\\jorda\\Desktop\\projects\\Electrum\\_data\\MIDAS-latest.txt')

    spec = importlib.util.spec_from_file_location("", os.path.abspath("./metaboverse_cli/analyze/model.py"))
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    gather_synonyms = model.gather_synonyms

    output_file = 'C:\\Users\\jorda\\Desktop\\HSA-latest.eldb'

def __main__(
        args_dict,
        network,
        data,
        columns,
        species_id,
        output_file):
    """Generate graph object for visualization
    """

    print('Preparing metadata...')
    # Generate output file name
    graph_name = name_graph(
        output_file=output_file,
        species_id=species_id)

    reverse_genes = {v:k for k,v in network['ensembl_synonyms'].items()}
    protein_dictionary = uniprot_ensembl_reference(
        uniprot_reference=network['uniprot_synonyms'],
        ensembl_reference=reverse_genes)
    progress_feed(args_dict, "model", 1)

    chebi_dictionary = build_chebi_reference(
        chebi=network['chebi_mapper'],
        uniprot=network['uniprot_metabolites'])

    name_reference = build_name_reference(
        ensembl=network['ensembl_synonyms'],
        uniprot=network['uniprot_synonyms'])

    # add any mapping IDs
    # Add synonyms
    # Change name to user provided if available
    metabolite_mapper = load_metabolite_synonym_dictionary()

    u = {}
    for k,v in network['uniprot_metabolites'].items():
        u[v] = k

    common_metabolites = get_metabolites(
        data=data,
        columns=columns)

    metaboanalyst_chebi_table = fetch_chebi_mappers(
        metabolites=common_metabolites)

    metabolites = define_mapper(
        mapper=metaboanalyst_chebi_table,
        data=data,
        columns=columns)

    # Generate graph
    # Name mapping
    print('Building network...')
    G = targeted_graph(
        metabolites=metabolites,
        reactions=network['reaction_database'],
        pathways=network['pathway_database'],
        species_reference=network['species_database'],
        name_database=network['name_database'],
        metabolite_mapper=metabolite_mapper,
        uniprot_mapper=u,
        component_database=network['components_database'])
    progress_feed(args_dict, "model", 9)

    degree_dictionary = compile_node_degrees(
        graph=G)

    output_database = json_graph.node_link_data(G)
    output_database['pathway_dictionary'] = network['pathway_database']
    output_database['degree_dictionary'] = degree_dictionary
    output_database['curation_date'] = date.today().strftime('%Y-%m-%d')

    with open(graph_name, 'w') as f:
        json.dump(output_database, f, indent=4)

    print('Graphing complete.')
    progress_feed(args_dict, "graph", 2)

    return graph_name