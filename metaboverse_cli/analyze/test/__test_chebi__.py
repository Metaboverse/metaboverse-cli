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
import zipfile
import importlib.util
import numpy as np
import networkx as nx
import xml.etree.ElementTree as et
import pandas as pd
import pickle
import os

print("Testing prepare_data.py")
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
prepare_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare_data)

zipped_net = os.path.abspath(
    './metaboverse_cli/analyze/test/HSA.zip')
with zipfile.ZipFile(zipped_net, 'r') as zip_file:
    zip_file.extractall(
        os.path.abspath(
            './metaboverse_cli/analyze/test'))

network_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/HSA.mvdb")
with open(network_url, 'rb') as network_file:
    network = pickle.load(network_file)

# CHEBI mapping
print('Testing analyze/__main__.py for CHEBI mapping...')
metabolomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/mixed_chebi_ids.txt")
__main__ = prepare_data.__main__
data, stats, unmapped = __main__(
    network=network,
    transcriptomics_url="None",
    proteomics_url="None",
    metabolomics_url=metabolomics_url)
assert data.shape == (8, 1), "Mixed CHEBI ID mapping failed"
"""
CHEBI:57972	                -0.187287     ->  "L-Ala"
CHEBI:18012	                -0.200046     ->  "FUMA"
CHEBI:30797	                -0.108502     ->  "MAL"
D-glucose	                -0.424363     ->  "Glc"
Fructose-6-phosphate	     0.066205     ->  "Fru(6)P"
L-LacTic Acid	            -0.000151     ->  "LACT"
SuCcINic acId	             0.068246     ->  "SUCCA"
gibberish	                -0.110443     ->  N/A
"""

args_chebi = {
    'database_source': 'reactome',
    'network': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA.mvdb'),
    'organism_curation': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA.mvdb'),
    'organism_id': 'HSA',
    'transcriptomics': 'None',
    'proteomics': 'none',
    'metabolomics': metabolomics_url,
    'output_file': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA_test.mvrs'),
    'collapse_with_modifiers': False,
    'broadcast_genes': True,
    'broadcast_metabolites': True,
    'labels': '0',
    'blocklist': ''
}
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/__main__.py"))
__main__ = importlib.util.module_from_spec(spec)
spec.loader.exec_module(__main__)
test_modeling = __main__.__main__
test_modeling(args_chebi)

with open(args_chebi['output_file']) as f:
    chebi_json = json.load(f)

for n in chebi_json['nodes']:
    if n['name'] == 'L-Ala':
        assert n['values'] == [-0.187286934], "Mixed CHEBI mapping failed"
    if n['name'] == 'SUCCA':
        assert n['values'] == [0.068245646], "Mixed CHEBI mapping failed"
    if n['values'] == [-0.11044283]:
        raise Exception("Mixed CHEBI mapping failed")
    if n['name'] == 'Fru(6)P':
        assert n['values'] == [0.06620505], "Mixed CHEBI mapping failed"
    if n['name'] == 'Glc':
        assert n['values'] == [-0.424362985], "Mixed CHEBI mapping failed"
    if n['name'] == 'LACT':
        assert n['values'] == [-0.00015059899999999999], "Mixed CHEBI mapping failed"
    if n['name'] == 'Glc':
        assert n['values'] == [-0.424362985], "Mixed CHEBI mapping failed"
    if n['name'] == 'FUMA':
        assert n['values'] == [-0.200045781], "Mixed CHEBI mapping failed"
    if n['name'] == 'MAL':
        assert n['values'] == [-0.10850223], "Mixed CHEBI mapping failed"

os.remove(args_chebi['output_file'])
os.remove(args_chebi['network'])

print('Tests completed')
