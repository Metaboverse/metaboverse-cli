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
import importlib.util
import pickle
import os

spec = importlib.util.spec_from_file_location(
    "__main__", os.path.abspath("./metaboverse_cli/mapper/__main__.py"))
mapper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mapper)

# Run
args_dict = {
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "mapper", "test")) + os.path.sep}
mapper.__main__(
    args_dict)

# Check
mapper_url = os.path.abspath(
    os.path.join(".", "metaboverse_cli", "mapper", "test", "metabolite_mapping.pickle"))
with open(mapper_url, 'rb') as mapper_file:
    mapper = pickle.load(mapper_file)

assert list(mapper.keys()) == ['hmdb_dictionary', 'display_dictionary',
                               'mapping_dictionary'], 'metabolite mapper failed to generate dictionaries'

key0 = list(mapper['hmdb_dictionary'].keys())[0]
assert type(mapper['hmdb_dictionary'][key0]
            ) == list, 'HMDB dictionary improperly formatted'

key1 = list(mapper['display_dictionary'].keys())[0]
assert type(mapper['display_dictionary'][key1]
            ) == list, 'Display dictionary improperly formatted'

key2 = list(mapper['mapping_dictionary'].keys())[0]
assert type(mapper['mapping_dictionary'][key2]
            ) == str, 'Mapping dictionary improperly formatted'

# Clean
os.remove(
    os.path.abspath(
        os.path.join(".", "metaboverse_cli", "mapper", "test", "metabolite_mapping.pickle")))

print('Tests completed')
