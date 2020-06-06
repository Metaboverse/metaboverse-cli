"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
    Copyright (C) 2019  Jordan A. Berg
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
import os
import json

"""Main/Utils
"""
# Functions to test
from metaboverse_cli.utils import update_session, \
                                    progress_feed, \
                                    check_directories, \
                                    check_files, \
                                    check_curate, \
                                    check_analyze, \
                                    argument_checks

# update_session()
session_file = os.path.abspath("./tests/test_data/session_data.json")
update_session(
    session_file=session_file,
    key='database_url',
    value='this is a fake')
with open(session_file) as json_file:
    data = json.load(json_file)
    assert data['database_url'] == 'this is a fake', 'update_session() failed'

update_session(
    session_file=session_file,
    key='database_url',
    value='')
with open(session_file) as json_file:
    data = json.load(json_file)
    assert data['database_url'] == '', 'update_session() failed'

# progress_feed()
progress_file = os.path.abspath("./tests/test_data/progress_data.json")
args_dict = {'progress_log': progress_file}
progress_feed(
    args_dict=args_dict,
    process="tester",
    amount=1)
with open(progress_file) as json_file:
    data = json.load(json_file)
    assert data['tester'] == 1, 'progress_feed() failed'

progress_feed(
    args_dict=args_dict,
    process="tester",
    amount=-1)
with open(progress_file) as json_file:
    data = json.load(json_file)
    assert data['tester'] == 0, 'progress_feed() failed'

# check_directories()
args_dict = {'output': os.path.abspath("./tests/test_data")}
dir = check_directories(
    input=args_dict['output'],
    argument='output')
assert dir == os.path.abspath("./tests/test_data") + '/', 'check_directories() failed'

# check_files()
args_dict = {'input': os.path.abspath("./tests/test_data/rnaseq_test.txt")}
f = check_files(
    input=args_dict['input'],
    argument='input')
assert f == os.path.abspath("./tests/test_data/rnaseq_test.txt"), 'check_files() failed'

# check_curate()
args_dict = {
    'output': os.path.abspath("./tests/test_data"),
    'species_id': 'SCE'}
try:
    check_curate(
        args_dict=args_dict)
except:
    raise Exception('check_curate() failed')

# argument_checks()
args_dict = {
    'output': os.path.abspath("./tests/test_data"),
    'species_id': 'SCE',
    'progress_log': progress_file,
    'session_data': session_file,
    'cmd': "BAD"}
args_dict = argument_checks(
    args_dict=args_dict)
assert args_dict == {
    'output': '/Users/jordan/Desktop/metaboverse-cli/tests/test_data/',
    'species_id': 'SCE',
    'progress_log': '/Users/jordan/Desktop/metaboverse-cli/tests/test_data/progress_data.json',
    'session_data': '/Users/jordan/Desktop/metaboverse-cli/tests/test_data/session_data.json',
    'cmd': 'BAD'}, 'argument_checks() failed'

"""Curation/Utils
"""
from metaboverse_cli import __main__ as curate


# test __main__() -- functional test
args_dict = {
    'species': 'SCE',
    'output': os.path.abspath("./tests/test_data")}

curate(
    args_dict=args_dict)

with open(args_dict['output'] + 'SCE_metaboverse_db.pickle', 'rb') as network_file:
    reactome_database = pickle.load(network_file)

assert reactome_database['master_reference']['R-ALL-389536'] == 'CO2', 'Unable to extract element from master reference'
assert reactome_database['pathways']['R-HSA-2562578']['reactions']['R-HSA-2562541']['name'] == 'TLR4-induced ripoptosome assembly', 'Unable to extract reaction name'

os.remove(args_dict['output'] + 'SCE_metaboverse_db.pickle')







"""Model/Utils
"""
