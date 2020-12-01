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
import os
import io
import pickle
import zipfile
import requests
import xml.etree.ElementTree as et
import pandas as pd

def write_database(
        output,
        file,
        database):
    """Write reactions database to pickle file
    """

    # Check provided path exists
    if not os.path.isdir(output):
        os.makedirs(output)

    # Clean up path
    dir = os.path.abspath(output) + '/'

    # Write information to file
    with open(dir + file, 'wb') as file_product:
        pickle.dump(database, file_product)

def parse_hmdb_synonyms(
        output_dir,
        url='https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
        file_name='hmdb_metabolites',
        xml_tag = '{http://www.hmdb.ca}'):
    """Retrieve HMDB chemical entity synonyms
    """

    output_file = output_dir + file_name
    print("Downloading HMDB metabolite reference...")
    hmdb_url = requests.get(url)
    if hmdb_url.ok:
        print("Unzipping HMDB metabolite reference...")
        hmdb_zip = zipfile.ZipFile(io.BytesIO(hmdb_url.content))
        hmdb_zip.extractall(output_dir)
        hmdb_zip = None
    else:
        raise Exception("Unable to download file at: " + url)

    print("Parsing HMDB metabolite records...")
    hmdb_contents = et.parse(output_file + '.xml')
    contents = hmdb_contents.getroot()
    os.remove(output_file + '.xml')

    hmdb_dictionary = {}
    display_dictionary = {}
    mapping_dictionary = {}

    for x in contents:

        if x.tag == xml_tag + 'metabolite':

            name = ''
            synonyms = set()
            display_synonyms = set()

            for y in x:

                if y.text != None:
                    simple_string = ''.join(
                        str(c).lower() for c in y.text if c.isalnum()
                    )
                else:
                    simple_string = 'None'

                if y.tag == xml_tag + 'name':
                    name = simple_string
                    synonyms.add(simple_string)
                    synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))

                if y.tag == xml_tag + 'synonyms':
                    for child in y:
                        simple_child = ''.join(
                            str(c).lower() for c in child.text if c.isalnum()
                        )
                        synonyms.add(simple_child.lower())
                        synonyms.add(str(child.text).lower())
                        display_synonyms.add(str(child.text))

                if y.tag == xml_tag + 'iupac_name':
                    synonyms.add(simple_string)
                    synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))

                if y.tag == xml_tag + 'traditional_iupac':
                    synonyms.add(simple_string)
                    synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))

            if name != '':
                hmdb_dictionary[name] = sorted(list(synonyms))
                display_dictionary[name] = sorted(list(display_synonyms))
                for l in hmdb_dictionary[name]:
                    mapping_dictionary[l] = name

    return hmdb_dictionary, display_dictionary, mapping_dictionary

def __main__(
        args_dict):
    """Build metabolite name mapping dictionary
    """

    hmdb_dictionary, display_dictionary, mapping_dictionary = parse_hmdb_synonyms(
        output_dir=args_dict['output']
    )

    mapping_db = {
        'hmdb_dictionary': hmdb_dictionary,
        'display_dictionary': display_dictionary,
        'mapping_dictionary': mapping_dictionary
    }

    print('Writing database to file...')
    write_database(
        output=args_dict['output'],
        file='metabolite_mapping.pickle',
        database=mapping_db
    )
