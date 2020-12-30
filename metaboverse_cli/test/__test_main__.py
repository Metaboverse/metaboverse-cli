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

"""Import dependencies
"""
import os
import pickle
import pandas as pd
import xml.etree.ElementTree as et
import networkx as nx
import numpy as np
import importlib.util

# Run full test on data
print("Testing analyze/__main__.py for modeling data")
spec = importlib.util.spec_from_file_location("", os.path.abspath("./metaboverse_cli/analyze/__main__.py"))
__main__ = importlib.util.module_from_spec(spec)
spec.loader.exec_module(__main__)
test_modeling = __main__.__main__

args_dict = {
    'database_source': 'reactome',
    'network': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA.mvdb'),
    'organism_curation': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA.mvdb'),
    'organism_id': 'HSA',
    'transcriptomics': os.path.abspath(
        './metaboverse_cli/analyze/test/rna_mapping_test.txt'),
    'proteomics': 'none',
    'metabolomics': os.path.abspath(
        './metaboverse_cli/analyze/test/metabolite_mapping_test.txt'),
    'output_file': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA_test.mvrs'),
    'collapse_with_modifiers': False,
    'broadcast_genes': True,
    'labels': '0',
    'blocklist': ''
}

test_modeling(args_dict)

rna_unmapped = os.path.abspath(
    './metaboverse_cli/analyze/test/rna_mapping_test_unmapped.txt'
)
rna = pd.read_csv(rna_unmapped, sep='\t', index_col=0)
assert len(rna.index.tolist()) == 7029, 'RNA mapping experienced error'
os.remove(rna_unmapped)

metabolite_unmapped = os.path.abspath(
    './metaboverse_cli/analyze/test/metabolite_mapping_test_unmapped.txt'
)
met = pd.read_csv(metabolite_unmapped, sep='\t', index_col=0)
assert met.index.tolist() == ['bMethyl.2.oxovalerate', 'DSS', 'Phenylacetylglycine'], 'Metabolite mapping experienced error'
os.remove(metabolite_unmapped)

os.remove(args_dict['output_file'])

print('Tests completed')
