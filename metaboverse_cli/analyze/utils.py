"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

Copyright (C) Jordan A. Berg

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
import pandas as pd
import os


def file_path(
        input):
    # Check input is contains full path address

    return os.path.abspath(input)


def check_suffix(
        file):
    # Get file suffix

    if file.split('.')[-1] == 'csv':
        suffix = ','
    elif file.split('.')[-1] == 'tsv':
        suffix = '\t'
    elif file.split('.')[-1] == 'txt':
        suffix = '\t'
    else:
        raise Exception(
            'Invalid data file provided. Expected a tab- or comma-delimited file')

    return suffix


def add_data(
        file):
    # Input data type
    # Check that file has full path
    file = file_path(file)

    # Figure out file type
    suffix = check_suffix(file)

    # Import dataframe
    data = pd.read_csv(
        file,
        sep=suffix,
        header=0,
        index_col=0,
        low_memory=False)

    return data


def convert_rgba(
        rgba_tuples,
        N=255):
    """Convert python RGBA tuple to web-friendly tuple for later viz
    """

    js = []
    for x in rgba_tuples:

        rgba_list = list(x)
        rgba_new = []
        for x in rgba_list[:3]:
            rgba_new.append(int(x * N))

        rgba_new.append(rgba_list[3])

        js.append(tuple(rgba_new))

    return js


def remove_defective_reactions(
        network):
    """
    """
    no_defective_reactions = {}
    for key in network['reaction_database'].keys():
        rxn_name = network['reaction_database'][key]['name'].lower()
        if 'defective' not in rxn_name \
                and 'mutant' not in rxn_name:
            no_defective_reactions[key] = network['reaction_database'][key]

    return no_defective_reactions
