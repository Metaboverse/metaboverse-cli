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
import sys
import re
import stat
import tarfile
import time
import hashlib
import xml.etree.ElementTree as et
import glob

"""Import internal dependencies
"""
try:
    from utils import progress_feed, update_session, safestr
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location("", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    update_session = utils.update_session
    safestr = utils.safestr

"""Global variables
"""
sbml_namespace = '{{http://www.sbml.org/sbml/level{0}/version{1}/core}}'
rdf_namespace = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
bqbiol_namespace = '{http://biomodels.net/biology-qualifiers/}'
sbml_level = '3'
sbml_version = '1'

"""Functions
"""
def test():

    output_dir = '/Users/jordan/Desktop/'
    pathways_dir = '/Users/jordan/Desktop/metaboverse_data/records/HSA/'
    species_id = 'HSA'
    args_dict = None

def handle_folder_contents(dir):

    if os.path.exists(dir):
        try:
            files = glob.glob(dir + "*")
            for f in files:
                os.remove(f)
        except:
            print('Unable to remove files from: ' + str(dir) + ' ... skipping...')

        try:
            os.rmdir(dir)
        except:
            print('Unable to remove directory named: ' + str(dir) + ' ... skipping...')

def unpack_pathways(
        output_dir,
        url='https://reactome.org/download/current/all_species.3.1.sbml.tgz'):
    """Load tarballed sbml reactome pathway files from reactome site
    """

    file = output_dir + url.split('/')[-1]
    pathways_dir = file[:-4] + os.path.sep
    handle_folder_contents(
        dir=pathways_dir)

    os.system('curl -L ' + url + ' -o \"' + file + '\"')
    os.makedirs(pathways_dir)

    tar = tarfile.open(file, "r:gz")
    tar.extractall(path=pathways_dir)
    tar.close()
    os.remove(file)

    return pathways_dir

def get_pathways(
        species_id,
        pathways_dir):
    """Get list of pathways to parse
    """

    # Check provided path exists
    if not os.path.isdir(pathways_dir):
        raise Exception(pathways_dir, 'does not exist')

    # Clean up path
    if os.path.abspath(pathways_dir).endswith(os.path.sep):
        dir = os.path.abspath(pathways_dir)
    else:
        dir = os.path.abspath(pathways_dir) + os.path.sep

    # Get list of files and their reaction name
    file_list = os.listdir(dir)
    pathways_list = [f for f in file_list if species_id in f]
    pathways_list = [f.split('.')[:-1][0] for f in pathways_list]

    return pathways_list

def get_database(
        pathways_dir,
        pathway_name,
        extension='.sbml'):
    """Import sbml reaction data
    """

    if not pathways_dir.endswith(os.path.sep):
        pathways_dir = pathways_dir + os.path.sep

    pathway_file = pathways_dir + pathway_name + extension
    pathway_contents = et.parse(pathway_file)
    contents = pathway_contents.getroot()

    return contents

def get_metadata(
        reaction,
        sbml_level,
        sbml_version,
        sbml_namespace=sbml_namespace):
    """Get basic metadata for a reaction
    """

    compartment = safestr(reaction.attrib['compartment'])
    id = safestr(reaction.attrib['id'])
    name = safestr(reaction.attrib['name'])

    reversible = safestr(reaction.attrib['reversible'])
    if reversible == 'false':
        if '<' in name and '>' in name:
            reversible = 'true'

    try:
        notes = safestr(reaction.findall(
            str(sbml_namespace + 'notes').format(
                sbml_level,
                sbml_version
            )
        )[0][0].text)
    except:
        notes = ''
        print('No notes available for', name)

    return compartment, id, name, reversible, notes

def add_reaction(
        pathway_database,
        reaction,
        pathway,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Add reactions to pathway
    """

    _id = reaction.attrib['id']
    pathway_database[pathway]['reactions'].add(_id)

    return pathway_database, _id

def add_reaction_components(
        type,
        reaction,
        sbml_namespace=sbml_namespace,
        sbml_level=sbml_level,
        sbml_version=sbml_version):
    """Add reaction components to reactions database
    For type, options are "listOfReactants", "listOfProducts", or
    "listOfModifiers"
    """

    # Collect modifiers for a given reaction by species ID
    component_list = reaction.findall(
        str(sbml_namespace + type).format(
            sbml_level,
            sbml_version
        )
    )

    if len(component_list) > 0:
        component_list = component_list[0]

    items = []
    for child in component_list:

        if 'modifier' in child.attrib['id']:

            type = None
            if 'catalyst' in child.attrib['id'] \
            or 'positive' in child.attrib['id']:
                type = 'catalyst'

            elif 'inhibitor' in child.attrib['id'] \
            or 'negative' in child.attrib['id']:
                type = 'inhibitor'

            else:
                type = 'other'

            items.append([child.attrib['species'], type])

        else:
            items.append(child.attrib['species'])

    return items

def add_reaction_components_manual(
        type,
        reaction,
        sbml_namespace=sbml_namespace,
        sbml_level=sbml_level,
        sbml_version=sbml_version):
    """Add reaction components to reactions database
    For type, options are "listOfReactants", "listOfProducts", or
    "listOfModifiers"
    """

    # Collect modifiers for a given reaction by species ID
    component_list = reaction.findall(
        str(sbml_namespace + type).format(
            sbml_level,
            sbml_version
        )
    )

    if len(component_list) > 0:
        component_list = component_list[0]

    items = []
    for child in component_list:

        if type == 'listOfModifiers':
            _type = 'modifier'
            items.append([child.attrib['species'], _type])

        else:
            items.append(child.attrib['species'])

    return items

def add_names(
        name_database,
        child,
        specie,
        search_string='is',
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Add names to dictionary to map species ID
    """

    for rank in child.iter(str(bqbiol_namespace + search_string)):
        for _rank in rank.iter(str(rdf_namespace + 'li')):

            item = _rank.attrib[str(rdf_namespace + 'resource')]
            _id = item.split('/')[-1]
            if 'chebi' in item.lower():
                _id = check_chebi(item=_id)
                _id = _id.split(' ')[0]
            name_database[_id] = specie
            name_database[specie] = specie

            # If element has parentheses, remove what's in between as
            # additional key
            if '(' in _id and ')' in _id:
                name_database = add_alternative_names(
                    name_database=name_database,
                    item=_id,
                    specie=specie)

    return name_database

def add_bigg_names(
        name_database,
        child,
        specie,
        search_string='notes',
        sbml_namespace=sbml_namespace,
        sbml_level=sbml_level,
        sbml_version=sbml_version):
    """Add names to dictionary to map species ID
    """

    for c in child:
        if c.tag == str(sbml_namespace + 'notes').format(
                sbml_level,
                sbml_version):
            for z in c[0]:
                _z = z.text
                if 'bigg' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie
                if 'biopath' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie
                if 'kegg' in _z.lower():
                    if ',' in _z:
                        for _z_ in _z.split(','):
                            name_database[_z_.replace(' ', '')] = specie
                    else:
                        name_database[_z.replace(' ', '')] = specie
                if 'metacyc' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie
                if 'mxnref' in _z.lower():
                    name_database[_z.split(':')[1].replace(' ', '')] = specie

    name_database[specie] = specie
    return name_database

def add_alternative_names(
        name_database,
        item,
        specie):
    """Add alternative names to name database for mapping
    """

    _remove = item[item.find('(') : item.find(')') + 1]
    mod_item = item.replace(_remove, '')
    name_database[mod_item] = specie

    return name_database

def check_chebi(
        item):
    """Some special formatting handling for CHEBI entries
    """

    item_parsed = item.lower().split('chebi:')[1]
    item_returned = 'CHEBI:' + item_parsed

    return item_returned

def add_species(
        species_database,
        name_database,
        compartment_database,
        components_database,
        pathway_record,
        sbml_namespace=sbml_namespace,
        sbml_level=sbml_level,
        sbml_version=sbml_version,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Add species records for pathway to database
    """

    species = pathway_record.findall(
        str(sbml_namespace + 'listOfSpecies').format(
            sbml_level,
            sbml_version
        )
    )[0]

    # Collect species information
    for child in species:

        # Initialize specie record and add common name
        specie = child.attrib['id']
        name = child.attrib['name']
        compartment = child.attrib['compartment']

        if '[' in name:
            name = name.split(' [')[0]

        species_database[specie] = name

        # Add compartment membership for specie
        compartment_database[specie] = compartment

        # Add names and ids to name dictionary
        name_database[name] = specie

        components_database[specie] = {
            'id': specie,
            'reactome_id': '',
            'name': name,
            'is': '',
            'isEncodedBy':'',
            'hasPart': [],
            'type': '',
            'compartment': compartment
        }

        for rank in child.iter(str(bqbiol_namespace + 'is')):
            for _rank in rank.iter(str(rdf_namespace + 'li')):
                item = _rank.attrib[str(rdf_namespace + 'resource')]
                if 'reactome' not in item.lower():
                    if 'chebi' in item.lower():
                        _id = item.split('chebiId=')[1]
                        components_database[specie]['is'] = _id
                        components_database[specie]['type'] = 'metabolite_component'
                    elif 'uniprot' in item.lower():
                        _id = item.split('/')[-1]
                        components_database[specie]['is'] = _id
                        components_database[specie]['type'] = 'protein_component'
                    elif 'mirbase' in item.lower():
                        _id = item.split('acc=')[1]
                        components_database[specie]['is'] = _id
                        components_database[specie]['type'] = 'mirna_component'
                    else:
                        components_database[specie]['type'] = 'other'
                else:
                    r_id = item.split('/')[-1]
                    components_database[specie]['reactome_id'] = r_id

        for rank in child.iter(str(bqbiol_namespace + 'hasPart')):
            for _rank in rank.iter(str(rdf_namespace + 'li')):
                item = _rank.attrib[str(rdf_namespace + 'resource')]
                if 'reactome' not in item:
                    components_database[specie]['type'] = 'complex_component'
                    if 'chebi' in item.lower():
                        _id = item.split('chebiId=')[1]
                        components_database[specie]['hasPart'].append(_id)
                    elif 'uniprot' in item.lower():
                        _id = item.split('/')[-1]
                        components_database[specie]['hasPart'].append(_id)
                    elif 'mirbase' in item.lower():
                        _id = item.split('acc=')[1]
                        components_database[specie]['hasPart'].append(_id)
                    else:
                        pass

        # Add source ID
        name_database = add_names(
            name_database=name_database,
            child=child,
            specie=specie,
            search_string='is',
            bqbiol_namespace=bqbiol_namespace,
            rdf_namespace=rdf_namespace)

    return (species_database, name_database, compartment_database,
        components_database)

def process_components(
        output_dir,
        pathways_dir,
        pathways_list,
        species_id,
        args_dict=None,
        sbml_namespace=sbml_namespace,
        bqbiol_namespace=bqbiol_namespace,
        rdf_namespace=rdf_namespace):
    """Process species-specific pathways
    """

    # Initialize databases
    pathway_database = {}
    reaction_database = {}
    species_database = {}
    name_database = {}
    compartment_database = {}
    compartment_dictionary = {}
    components_database = {}

    print('Extracting pathway-level reaction data for: ' + str(species_id))

    # Cycle through each pathway database and extract  contents
    for pathway in pathways_list:

        db = get_database(
            pathways_dir,
            pathway)
        sbml_level = db.attrib['level']
        sbml_version = db.attrib['version']

        pathway_record = db.findall(
            str(sbml_namespace + 'model').format(
                sbml_level,
                sbml_version
            )
        )[0]

        pathway_info = pathway_record.attrib

        id = pathway_info['id']
        pathway_database[pathway] = {
            'id': id,
            'reactome': pathway,
            'name': pathway_info['name'],
            'reactions': set()
        }

        # Parse out reactions
        reactions = pathway_record.findall(
            str(sbml_namespace + 'listOfReactions').format(
                sbml_level,
                sbml_version
            )
        )[0]

        # Parse out compartment IDs and names
        compartments = pathway_record.findall(
            str(sbml_namespace + 'listOfCompartments').format(
                sbml_level,
                sbml_version
            )
        )[0]
        for c in range(len(compartments)):
            id = compartments[c].attrib['id']
            name = compartments[c].attrib['name']
            compartment_dictionary[id] = name

        # Extract reactions from pathway
        for reaction in reactions:

            # Get metadata
            compartment, id, name, reversible, notes = get_metadata(
                reaction=reaction,
                sbml_level=sbml_level,
                sbml_version=sbml_version,
                sbml_namespace=sbml_namespace)

            # Get pathway high-level information (reactions, name, compartment)
            pathway_database, reaction_id = add_reaction(
                pathway_database=pathway_database,
                reaction=reaction,
                pathway=pathway,
                bqbiol_namespace=bqbiol_namespace,
                rdf_namespace=rdf_namespace)

            name_database[name] = reaction_id
            reaction_database[id] = {
                'compartment': compartment,
                'id': id,
                'name': name,
                'reversible': reversible,
                'notes': notes}

            # Collect reactants for a given reaction by species ID
            reaction_database[reaction_id]['reactants'] = add_reaction_components(
                type='listOfReactants',
                reaction=reaction,
                sbml_namespace=sbml_namespace,
                sbml_level=sbml_level,
                sbml_version=sbml_version)

            # Collect products for a given reaction by species ID
            reaction_database[reaction_id]['products'] = add_reaction_components(
                type='listOfProducts',
                reaction=reaction,
                sbml_namespace=sbml_namespace,
                sbml_level=sbml_level,
                sbml_version=sbml_version)

            # Collect modifiers for a given reaction by species ID
            reaction_database[reaction_id]['modifiers'] = add_reaction_components(
                type='listOfModifiers',
                reaction=reaction,
                sbml_namespace=sbml_namespace,
                sbml_level=sbml_level,
                sbml_version=sbml_version)

        # Convert reaction set for pathway to list
        pathway_database[pathway]['reactions'] = list(
            pathway_database[pathway]['reactions'])

        # Generate species dict
        species_database, name_database, compartment_database, \
        components_database = add_species(
            species_database=species_database,
            name_database=name_database,
            compartment_database=compartment_database,
            components_database=components_database,
            pathway_record=pathway_record,
            sbml_namespace=sbml_namespace,
            sbml_level=sbml_level,
            sbml_version=sbml_version,
            bqbiol_namespace=bqbiol_namespace,
            rdf_namespace=rdf_namespace)

    return (pathway_database, reaction_database, species_database,
        name_database, compartment_database, compartment_dictionary,
        components_database)

def load_sbml(
        sbml_url):
    """Load supported SBML file for organism network curation
    """

    parsed_url = sbml_url.split(os.path.sep)
    parsed_path = os.path.sep.join(parsed_url[0:-1]) + os.path.sep
    parsed_file = parsed_url[-1].split('.')[0]
    parsed_extension = '.' + parsed_url[-1].split('.')[1]

    return get_database(
        parsed_path,
        parsed_file,
        extension=parsed_extension)

def update_model_metadata(
        sbml_db,
        args_dict):
    """Get model metadata and update session info
    """

    session_file = args_dict['session_data']
    update_session(
        session_file=session_file,
        key='organism_id',
        value=sbml_db[0].attrib['id'])
    if 'name' in sbml_db[0].attrib:
        update_session(
            session_file=session_file,
            key='organism',
            value=sbml_db[0].attrib['name'])
    else:
        update_session(
            session_file=session_file,
            key='organism',
            value='unknown')
    if 'metaid' in sbml_db[0].attrib:
        update_session(
            session_file=session_file,
            key='database_version',
            value=sbml_db[0].attrib['metaid'] + ' (' + args_dict['database_source'] + ')')
    else:
        update_session(
            session_file=session_file,
            key='database_version',
            value='unknown')

def process_manual(
        sbml_db,
        args_dict=None):
    """Parse network curation elements
    """

    if 'core' in sbml_db.tag:
        _a = '/core'
    else:
        _a = ''

    sbml_namespace = '{{http://www.sbml.org/sbml/level{0}/version{1}' + _a + '}}'
    sbml_level = sbml_db.attrib['level']
    sbml_version = sbml_db.attrib['version']

    # Initialize databases
    pathway_database = {
        'All': {
            'id': 'All',
            'reactome': 'All',
            'name': 'All',
            'reactions': set()
        }
    }
    reaction_database = {}
    name_database = {}
    compartment_dictionary = {}
    compartment_database = {}
    species_database = {}
    components_database = {}

    # Get model information
    if args_dict != None:
        update_model_metadata(
            sbml_db=sbml_db,
            args_dict= args_dict
        )

    # Get model categories
    elements = [x for x in sbml_db[0]]

    # Generate compartment dictionary
    for x in elements:
        if x.tag == str(sbml_namespace + 'listOfCompartments').format(
                sbml_level,
                sbml_version):
            for child in x:
                id = child.attrib['id']
                name = child.attrib['name']
                compartment_dictionary[id] = name

    #Generate species database
    for x in elements:
        if x.tag == str(sbml_namespace + 'listOfSpecies').format(
                sbml_level,
                sbml_version):
            for child in x:
                specie = child.attrib['id']
                if 'name' in child.attrib:
                    name = child.attrib['name']
                else:
                    name = specie
                if 'sboTerm' in child.attrib:
                    sboTerm = child.attrib['sboTerm']
                else:
                    sboTerm = ''
                compartment = child.attrib['compartment']

                species_database[specie] = name
                compartment_database[specie] = compartment
                name_database[name] = specie
                components_database[specie] = {
                    'id': specie,
                    'reactome_id': sboTerm,
                    'name': name,
                    'is': specie,
                    'isEncodedBy':'',
                    'hasPart': [],
                    'type': '',
                    'compartment': compartment
                }

                for rank in child.iter(str(bqbiol_namespace + 'is')):
                    for _rank in rank.iter(str(rdf_namespace + 'li')):
                        item = _rank.attrib[str(rdf_namespace + 'resource')]
                        if 'reactome' not in item.lower():
                            if 'chebi' in item.lower() \
                            or 'kegg' in item.lower() \
                            or 'hmdb' in item.lower() \
                            or 'bigg' in item.lower():
                                _id = item.split('/')[-1]
                                components_database[specie]['is'] = _id
                                components_database[specie]['type'] = 'metabolite_component'
                            elif 'uniprot' in item.lower():
                                _id = item.split('/')[-1]
                                components_database[specie]['is'] = _id
                                components_database[specie]['type'] = 'protein_component'
                            else:
                                components_database[specie]['type'] = 'other'
                        else:
                            r_id = item.split('/')[-1]
                            components_database[specie]['reactome_id'] = r_id

                for rank in child.iter(str(bqbiol_namespace + 'hasPart')):
                    for _rank in rank.iter(str(rdf_namespace + 'li')):
                        item = _rank.attrib[str(rdf_namespace + 'resource')]
                        if 'reactome' not in item:
                            components_database[specie]['type'] = 'complex_component'
                            if 'chebi' in item.lower() \
                            or 'kegg' in item.lower() \
                            or 'hmdb' in item.lower() \
                            or 'bigg' in item.lower():
                                _id = item.split('/')[-1]
                                components_database[specie]['hasPart'].append(_id)
                            elif 'uniprot' in item.lower():
                                _id = item.split('/')[-1]
                                components_database[specie]['hasPart'].append(_id)
                            elif 'mirbase' in item.lower():
                                _id = item.split('acc=')[1]
                                components_database[specie]['hasPart'].append(_id)
                            else:
                                pass

                for rank in child.iter(str(bqbiol_namespace + 'isEncodedBy')):
                    for _rank in rank.iter(str(rdf_namespace + 'li')):
                        item = _rank.attrib[str(rdf_namespace + 'resource')]
                        if 'reactome' not in item:
                            if 'kegg.genes' in item.lower():
                                _id = item.split('/')[-1]
                                _id_ = _id.split(':')[-1]
                                components_database[specie]['isEncodedBy'] = _id_
                            else:
                                _id = item.split('/')[-1]
                                components_database[specie]['isEncodedBy'] = _id

                # Add source ID
                name_database = add_names(
                    name_database=name_database,
                    child=child,
                    specie=specie,
                    search_string='is',
                    bqbiol_namespace=bqbiol_namespace,
                    rdf_namespace=rdf_namespace)

                name_database = add_bigg_names(
                    name_database=name_database,
                    child=child,
                    specie=specie,
                    search_string='notes',
                    sbml_namespace=sbml_namespace,
                    sbml_level=sbml_level,
                    sbml_version=sbml_version)

    #Generate reaction database
    for x in elements:
        if x.tag == str(sbml_namespace + 'listOfReactions').format(
                sbml_level,
                sbml_version):
            for child in x:
                # Get metadata
                _id = child.attrib['id']
                if 'name' in child.attrib:
                    _name = child.attrib['name']
                else:
                    _name = _id
                if 'reversible' in child.attrib:
                    _reversible = child.attrib['reversible']
                else:
                    _reversible = 'false'

                name_database[_name] = _id
                pathway_database['All']['reactions'].add(_id)
                reaction_database[_id] = {
                    'compartment': '',
                    'id': _id,
                    'name': _name,
                    'reversible': _reversible,
                    'notes': ''}
                reaction_database[_id]['reactants'] = add_reaction_components_manual(
                    type='listOfReactants',
                    reaction=child,
                    sbml_namespace=sbml_namespace,
                    sbml_level=sbml_level,
                    sbml_version=sbml_version)
                reaction_database[_id]['products'] = add_reaction_components_manual(
                    type='listOfProducts',
                    reaction=child,
                    sbml_namespace=sbml_namespace,
                    sbml_level=sbml_level,
                    sbml_version=sbml_version)
                reaction_database[_id]['modifiers'] = add_reaction_components_manual(
                    type='listOfModifiers',
                    reaction=child,
                    sbml_namespace=sbml_namespace,
                    sbml_level=sbml_level,
                    sbml_version=sbml_version)

    return (pathway_database, reaction_database, species_database,
        name_database, compartment_database, compartment_dictionary,
        components_database)

def __main__(
        species_id,
        output_dir,
        database_source,
        sbml_url,
        args_dict):
    """Fetch all reactions for a given organism
    """

    # Get pathways files
    if database_source.lower() == 'reactome':
        pathways_dir = unpack_pathways(
            output_dir=output_dir)
        progress_feed(args_dict, "curate", 10)

        pathways_list = get_pathways(
            species_id=species_id,
            pathways_dir=pathways_dir)
        progress_feed(args_dict, "curate", 7)

        # Get list of reaction files to use for populating database
        pathway_database, reaction_database, species_database, \
        name_database, compartment_database, compartment_dictionary, \
        components_database = process_components(
            output_dir=output_dir,
            pathways_dir=pathways_dir,
            pathways_list=pathways_list,
            species_id=species_id,
            args_dict=args_dict)
        progress_feed(args_dict, "curate", 5)

        if 'sbml' in pathways_dir:
            handle_folder_contents(
                dir=pathways_dir)
        else:
            print('Could not find SMBL file directory, skipping removal of this directory...')

    elif database_source.lower() == 'biomodels/bigg' and sbml_url != "None":
        sbml_db = load_sbml(
            sbml_url=sbml_url)
        progress_feed(args_dict, "curate", 10)

        pathway_database, reaction_database, species_database, \
        name_database, compartment_database, compartment_dictionary, \
        components_database = process_manual(
            sbml_db=sbml_db,
            args_dict=args_dict)
        progress_feed(args_dict, "curate", 13)

    else:
        raise Exception('Input database type not supported by Metaboverse. If you would like the database type included, please submit an issue at <https://github.com/Metaboverse/Metaboverse/issues>.')

    return (pathway_database, reaction_database, species_database,
        name_database, compartment_dictionary, components_database)
