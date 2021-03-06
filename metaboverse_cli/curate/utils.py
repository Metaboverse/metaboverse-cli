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
import pandas as pd
import sys
import os


def get_table(
        output_dir,
        url,
        column_names,
        organism='Homo sapiens',
        organism_key='organism'):
    """Get reactome table from web
    """

    # chebi_reactome_reactions
    file = unpack_table(
        url=url,
        output_dir=output_dir)

    if isinstance(column_names, list):
        header_type = None
    else:
        header_type = column_names

    data = pd.read_csv(
        file,
        sep='\t',
        header=header_type,
        low_memory=False)

    if isinstance(column_names, list) \
            or organism == None:
        data.columns = column_names
        data_organism = data.loc[data[organism_key] == organism]

    else:
        data_organism = data

    return data_organism


"""Open reactome table from web
"""


def unpack_table(
        url,
        output_dir='./'):

    file = output_dir + url.split('/')[-1]
    os.system('curl -L ' + url + ' -o "' + file + '"')

    return file
