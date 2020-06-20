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
import re
import pandas as pd
from math import floor

# Check input is contains full path address
def file_path(
        input):

    return os.path.abspath(input)

# Get file suffix
def check_suffix(
        file):

    if file.split('.')[-1] == 'csv':
        suffix = ','
    elif file.split('.')[-1] == 'tsv':
        suffix = '\t'
    elif file.split('.')[-1] == 'txt':
        suffix = '\t'
    else:
        raise Exception('Invalid data file provided. Expected a tab- or comma-delimited file')

    return suffix

# Input data type
def add_data(
        file):

    # Check that file has full path
    file = file_path(file)

    # Figure out file type
    suffix = check_suffix(file)

    # Import dataframe
    data = pd.read_csv(
        file,
        sep = suffix,
        header = 0,
        index_col = 0,
        low_memory = False)

    return data

def convert_rgba(
        rgba_tuples):
    """Convert python RGBA tuple to web-friendly tuple for later viz
    """

    js = []
    for x in rgba_tuples:

        rgba_list = list(x)
        rgba_new = []
        for x in rgba_list[:3]:
            rgba_new.append(int(x * 255))

        rgba_new.append(rgba_list[3])

        js.append(tuple(rgba_new))

    return js
