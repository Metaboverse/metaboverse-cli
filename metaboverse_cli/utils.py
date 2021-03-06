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
import json
import pickle
import sys
import os

def read_network(
        network_url):
    """Read in network from previous curation module
    - was provided as a URL to the file and saved to args_dict['network'] in
    "curate" sub-module
    """

    with open(network_url, 'rb') as network_file:
        network = pickle.load(network_file)

    return network


def prepare_output(
        output):
    """Get output directory prepared
    """

    # Check provided path exists
    if not os.path.isdir(output):
        os.makedirs(output)

    # Clean up path
    if os.path.abspath(output).endswith(os.path.sep):
        dir = os.path.abspath(output)
    else:
        dir = os.path.abspath(output) + os.path.sep

    return dir


def write_database(
        output,
        file,
        database):
    """Write reactions database to pickle file
    """

    dir = prepare_output(
        output=output)

    # Write information to file
    with open(dir + file, 'wb') as file_product:
        pickle.dump(database, file_product)


def write_database_json(
        output,
        file,
        database):
    """Write reactions database to JSON file
    """

    dir = prepare_output(
        output=output)

    # Write information to file
    with open(dir + file, 'w') as file_product:
        json.dump(database, file_product, indent=4)


def safestr(obj):
    """Covert ascii text if needed
    """

    return str(obj).encode('ascii', 'ignore').decode('ascii')


"""Update session information
"""


def update_session(
        session_file,
        key,
        value):

    if os.path.exists(str(session_file)) and str(session_file) != 'None':

        with open(session_file) as json_file:
            session = json.load(json_file)
            session[key] = value

        with open(session_file, 'w') as outfile:
            json.dump(session, outfile)

    else:
        print("Session file not found: " + str(session_file))


def get_session_value(
        session_file,
        key):

    if os.path.exists(str(session_file)) and str(session_file) != 'None':

        with open(session_file) as json_file:
            session = json.load(json_file)
            if key in session:
                return session[key]
            else:
                print("Cannot find specified key: " + key)
                return 'unknown'

    else:
        print("File at " + str(session_file) + " does not exist.")
        return 'unknown'


"""JS progress feed
"""


def progress_feed(
        args_dict=None,
        process=None,
        amount=1):

    if args_dict != None:
        if 'progress_log' in args_dict \
                and str(args_dict['progress_log']) != 'None':
            feed_file = args_dict['progress_log']

            if os.path.exists(feed_file) and process != None:

                with open(feed_file) as json_file:
                    data = json.load(json_file)
                    data[process] += amount

                with open(feed_file, 'w') as outfile:
                    json.dump(data, outfile)
    else:
        print('Could not access local variables during progress_feed() update.')


"""Check directory formatting
"""


def check_directories(
        input,
        argument):

    # Check that a file wasn't passed in
    if os.path.isdir(input) != True:
        raise Exception(safestr(argument) + ': ' +
                        safestr(input) + ' is not a directory')

    # Check input directory name is formatted correctly and fix if necessary
    input = safestr(os.path.abspath(input))

    if not input.endswith(os.path.sep):
        input += os.path.sep

    return input


"""Check file formatting
"""


def check_files(
        input,
        argument):

    # Check that a file wasn't passed in
    if os.path.isfile(input) != True:
        raise Exception(safestr(argument) + ': ' +
                        safestr(input) + ' is not a file')

    # Check input directory name is formatted correctly and fix if necessary
    input = safestr(os.path.abspath(input))

    return input


"""Check curation arguments
"""


def check_curate(
        args_dict):

    should_exit = False

    if args_dict['organism_id'] == None \
            or args_dict['organism_id'].lower() == 'none' \
            or args_dict['organism_id'].lower() == 'null':

        print('\nIncorrect species identifier provided: ' +
              safestr(args_dict['organism_id']))
        print('Please refer to https://reactome.org/ContentService/data/species/all for a valid list of organisms')
        should_exit = True

    if args_dict['output'] == None \
            or not os.path.exists(os.path.dirname(args_dict['output'])):

        print('\nIncorrect output parameter provided: ' +
              safestr(args_dict['output']))
        should_exit = True

    if should_exit == True:
        sys.exit(1)


"""Run general checks on arguments
Not sub-module-specific
"""


def argument_checks(
        args_dict):

    # Check output file
    if 'output' in args_dict \
            and args_dict['output'] == None:
        args_dict['output'] = os.getcwd()
    args_dict['output'] = safestr(args_dict['output'])

    if not args_dict['output'].endswith(os.path.sep):
        args_dict['output'] = args_dict['output'] + os.path.sep

    # Check user-provided directory formatting
    for key, value in args_dict.items():

        if key == 'cmd':
            pass

        elif os.path.isdir(str(value)) == True:
            args_dict[key] = check_directories(
                args_dict[key],
                key)

        elif os.path.isfile(str(value)) == True:
            args_dict[key] = check_files(
                args_dict[key],
                key)

        else:
            pass

    return args_dict
