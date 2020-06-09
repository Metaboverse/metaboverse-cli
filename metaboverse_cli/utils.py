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

Progress bar:
    From https://gist.github.com/vladignatyev/06860ec2040cb497f0f3
    The MIT License (MIT)
    Copyright (c) 2016 Vladimir Ignatev

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the Software
    is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
    PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
    FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
    OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
    OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import print_function

"""Import dependencies
"""
import os
import sys
import gc
import time
import datetime
from multiprocessing import cpu_count
from multiprocessing import Pool
import itertools
import json

"""Update session information
"""
def update_session(
        session_file,
        key,
        value):

    if os.path.exists(session_file) and session_file != None:

        with open(session_file) as json_file:
            session = json.load(json_file)
            session[key] = value

        with open(session_file, 'w') as outfile:
            json.dump(session, outfile)

    else:
        print("File at " + str(session_file) + " does not exist.")

"""JS progress feed
"""
def progress_feed(
        args_dict=None,
        process=None,
        amount=1):

    if args_dict != None:
        if 'progress_log' in args_dict and args_dict['progress_log'] != None:
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
        raise Exception(str(argument) + ': ' + str(input) + ' is not a directory')

    # Check input directory name is formatted correctly and fix if necessary
    input = os.path.abspath(input)

    if input.endswith('/'):
        pass
    else:
        input += '/'

    return input

"""Check file formatting
"""
def check_files(
        input,
        argument):

    # Check that a file wasn't passed in
    if os.path.isfile(input) != True:
        raise Exception(str(argument) + ': ' + str(input) + ' is not a file')

    # Check input directory name is formatted correctly and fix if necessary
    input = os.path.abspath(input)

    return input

"""Check curation arguments
"""
def check_curate(
        args_dict):

    should_exit = False

    if args_dict['species_id'] == None \
    or args_dict['species_id'].lower() == 'none' \
    or args_dict['species_id'].lower() == 'null':

        print('\nIncorrect species identifier provided: ' + args_dict['species_id'])
        print('Please refer to https://reactome.org/ContentService/data/species/all for a valid list of organisms')
        should_exit = True

    if args_dict['output'] == None \
    or not os.path.exists(os.path.dirname(args_dict['output'])):

        print('\nIncorrect output parameter provided: ' + args_dict['output'])
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
        args_dict['output'] = os.getcwd() + '/'

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
