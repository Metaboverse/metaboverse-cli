"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
    Copyright (C) 2019-2020  Jordan A. Berg
    jordan <dot> berg <at> biochem <dot> utah <dot> edu

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
import pickle

try:
    from metaboverse_cli.mapper.__main__ import __main__ as mapper
except:
    import importlib
    importlib.import_module('metaboverse_cli.mapper.__main__')



# Run
args_dict = {'output': os.path.abspath("./metaboverse_cli/mapper/test") + '/'}
mapper(
    args_dict)

# Check
mapper_url = os.path.abspath("./metaboverse_cli/mapper/test") + '/metabolite_mapping.pickle'
with open(mapper_url, 'rb') as mapper_file:
    mapper = pickle.load(mapper_file)

assert list(mapper.keys()) == ['hmdb_dictionary', 'display_dictionary', 'mapping_dictionary'], 'metabolite mapper failed to generate dictionaries'

key0 = list(mapper['hmdb_dictionary'].keys())[0]
assert type(mapper['hmdb_dictionary'][key0]) == list, 'HMDB dictionary improperly formatted'

key1 = list(mapper['display_dictionary'].keys())[0]
assert type(mapper['display_dictionary'][key1]) == list, 'Display dictionary improperly formatted'

key2 = list(mapper['mapping_dictionary'].keys())[0]
assert type(mapper['mapping_dictionary'][key2]) == str, 'Mapping dictionary improperly formatted'

# Clean
os.system(
    'rm '
    + os.path.abspath("./metaboverse_cli/mapper/test") + '/metabolite_mapping.pickle')
