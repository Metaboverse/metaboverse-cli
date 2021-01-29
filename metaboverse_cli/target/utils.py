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

COLUMN_REFERENCE = {
    'Metabolite': 0,
    'Query_protein': 1,
    'Protein_name': 2,
    'Gene_ID': 3,
    'Uniprot_ID': 4,
    'Protein_complex': 5,
    'log2_abundance': 6,
    'log2_abundance_corrected': 7,
    'met_mean': 8,
    'met_sd': 9,
    'p_value': 10,
    'q_value': 11,
    'Common_metabolite_name': 12,
    'MIDAS_ID': 13,
    'Pool': 14,
    'Screened_concentration_µM': 15,
    'KEGG_ID': 16,
    'HMDB_ID': 17,
    'SMILES_metabolite': 18,
    'KEGG_pathway_association': 19
}


def import_midas(
        filename,
        delimiter='\t'):
    """Import MIDAS database from flat format file
    """

    data = pd.read_csv(
        filename,
        sep=delimiter,
        encoding='unicode_escape')

    if len(data.columns.tolist()) != len(list(COLUMN_REFERENCE.keys())):
        raise Exception('Improperly formatted datatable provided: ', filename)

    return data, COLUMN_REFERENCE
