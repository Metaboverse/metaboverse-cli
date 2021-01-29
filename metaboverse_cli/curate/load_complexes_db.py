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
import os

"""Import internal dependencies
"""
try:
    from curate.utils import get_table
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "get_table", os.path.abspath("./metaboverse_cli/curate/utils.py"))
    get_table = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(get_table)
    get_table = get_table.get_table

"""Get tables
"""


def __main__(
        output_dir):

    complex_participants = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/ComplexParticipantsPubMedIdentifiers_human.txt',
        column_names=0)
    os.remove(output_dir + 'ComplexParticipantsPubMedIdentifiers_human.txt')

    complex_pathway = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/Complex_2_Pathway_human.txt',
        column_names=0)
    os.remove(output_dir + 'Complex_2_Pathway_human.txt')

    return {
        'complex_participants': complex_participants,
        'complex_pathway': complex_pathway}
