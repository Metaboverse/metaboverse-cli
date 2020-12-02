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
import os
import sys
import json

"""Main/Utils
"""
# Functions to test
try:
    from metaboverse_cli.utils import update_session, \
                                        progress_feed, \
                                        check_directories, \
                                        check_files, \
                                        check_curate, \
                                        argument_checks, \
                                        get_session_value
except:
    from utils import update_session, \
                        progress_feed, \
                        check_directories, \
                        check_files, \
                        check_curate, \
                        argument_checks, \
                        get_session_value

# update_session()
session_file = os.path.abspath(os.path.join(".", "metaboverse_cli", "test", "session_data.json"))
update_session(
    session_file=session_file,
    key='database_url',
    value='this is a fake')
with open(session_file) as json_file:
    data = json.load(json_file)
    assert data['database_url'] == 'this is a fake', 'update_session() failed'

update_session(
    session_file="bad_path",
    key='database_url',
    value='this is a fake')

update_session(
    session_file=session_file,
    key='database_url',
    value='')
with open(session_file) as json_file:
    data = json.load(json_file)
    assert data['database_url'] == '', 'update_session() failed'

# progress_feed()
progress_file = os.path.abspath("./metaboverse_cli/test/progress_data.json")
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
args_dict = {
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "test"))}
dir = check_directories(
    input=args_dict['output'],
    argument='output')
assert dir == os.path.abspath(
    os.path.join(".", "metaboverse_cli", "test")
) + os.path.sep, 'check_directories() failed'

try:
    dir = check_directories(
        input="bad_path",
        argument='output')
except:
    pass
else:
    raise Exception('check_directories() failed')

# check_files()
args_dict = {
    'input': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "test", "rnaseq_test.txt"))}
f = check_files(
    input=args_dict['input'],
    argument='input')
assert f == os.path.abspath(
    os.path.join(".", "metaboverse_cli", "test", "rnaseq_test.txt")
), 'check_files() failed'

try:
    f = check_files(
        input="bad_file",
        argument='input')
except:
    pass
else:
    raise Exception('check_files() failed')

# check_curate()
args_dict = {
    'output': os.path.abspath(os.path.join(".", "metaboverse_cli", "test")),
    'organism_id': 'SCE'}
try:
    check_curate(
        args_dict=args_dict)
except:
    raise Exception('check_curate() failed')

args_dict = {
    'output': os.path.abspath(os.path.join(".", "metaboverse_cli", "test")),
    'organism_id': 'None'}
try:
    check_curate(
        args_dict=args_dict)
except:
    pass
else:
    raise Exception('check_curate() failed')

args_dict = {
    'output': "bad_path",
    'organism_id': 'SCE'}
try:
    check_curate(
        args_dict=args_dict)
except:
    pass
else:
    raise Exception('check_curate() failed')

# argument_checks()
args_dict = {
    'output': os.path.abspath(os.path.join(".", "metaboverse_cli", "test")),
    'organism_id': 'SCE',
    'progress_log': progress_file,
    'session_data': session_file,
    'cmd': "BAD"}
args_dict = argument_checks(
    args_dict=args_dict)
assert args_dict == {
    'output': os.path.abspath(os.path.join(".", "metaboverse_cli", "test")) + os.path.sep,
    'organism_id': 'SCE',
    'progress_log': os.path.abspath(os.path.join(".", "metaboverse_cli", "test")) + os.path.sep + 'progress_data.json',
    'session_data': os.path.abspath(os.path.join(".", "metaboverse_cli", "test")) + os.path.sep + 'session_data.json',
    'cmd': 'BAD'}, 'argument_checks() failed'

args_dict['output'] = None
args_dict = argument_checks(
    args_dict=args_dict)
assert args_dict['output'] != None, 'argument_checks() failed'

# get_session_value()
val1 = get_session_value(
    session_file=args_dict['session_data'],
    key="database_url")
assert val1 == "", 'get_session_value() failed'

val2 = get_session_value(
    session_file=args_dict['session_data'],
    key="database_badURL")
assert val2 == "unknown", 'get_session_value() failed'

val3 = get_session_value(
    session_file= os.path.abspath(os.path.join(".", "metaboverse_cli", "test")) + os.path.sep + 'bad_file.json',
    key="database_url")
assert val3 == "unknown", 'get_session_value() failed'

print('Tests completed')
