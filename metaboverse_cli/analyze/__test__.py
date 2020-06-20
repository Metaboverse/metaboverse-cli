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
import pickle
import pandas as pd
import xml.etree.ElementTree as et

"""prepare_data.py
"""
import importlib.util
spec = importlib.util.spec_from_file_location("", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
prepare_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare_data)
read_data = prepare_data.read_data
format_data = prepare_data.format_data
format_metabolomics = prepare_data.format_metabolomics
output_unmapped = prepare_data.output_unmapped
extract_data = prepare_data.extract_data
broadcast_transcriptomics = prepare_data.broadcast_transcriptomics
copy_columns = prepare_data.copy_columns
catenate_data = prepare_data.catenate_data
__main__ = prepare_data.__main__

network_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/HSA_metaboverse_db.pickle")
with open(network_url, 'rb') as network_file:
    network = pickle.load(network_file)

transcriptomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/transcriptomics.txt")
proteomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/proteomics.txt")
metabolomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/metabolomics.txt")

# read_data()
transcriptomics_df = read_data(
    url=transcriptomics_url)
assert len(transcriptomics_df.columns.tolist()) == 2, "read_data() failed"

# format_data()
e_sym = {}
for k, v in network['ensembl_synonyms'].items():
    e_sym[v] = k
    e_sym[k] = k
t_mapped, t_unmapped = format_data(
    data=transcriptomics_df,
    reference=e_sym)
assert t_mapped.shape[0] >= t_unmapped.shape[0], "format_data()"
assert t_mapped.shape[1] == 2, "format_data()"
assert t_unmapped.shape[1] == 2, "format_data()"
assert transcriptomics_df.loc['GCN1'][0] == t_mapped.loc['ENST00000300648'][0], 'format_data() failed'
try:
    t_mapped.loc['GCN1']
except KeyError:
    pass
else:
    raise Exception('format_data() failed')

# output_unmapped()
output_unmapped(
    data=t_unmapped,
    url=transcriptomics_url)

# extract_data()
_v, _s = extract_data(
    data=t_mapped)
assert _v.shape[1] == 1, "extract_data()"
assert _s.shape[1] == 1, "extract_data()"

# broadcast_transcriptomics()
# How is this mapped gene vs protein in graphing?
_p, _p_stats = broadcast_transcriptomics(
    transcriptomics=_v,
    transcriptomics_stats=_s,
    gene_dictionary=network['ensembl_synonyms'],
    protein_dictionary=network['uniprot_synonyms'])
assert _v.loc['ENST00000300648'][0] == _p.loc['Q92616'][0], 'broadcast_transcriptomics() failed'

# copy_columns()
data_col, stat_col = copy_columns(_v, _s, 6)
assert len(_v.columns.tolist()) == 1, 'copy_columns() failed'
assert len(data_col.columns.tolist()) == 6, 'copy_columns() failed'
assert data_col[0].values.all() == data_col[5].values.all(), 'copy_columns() failed'
assert stat_col[0].values.all() == stat_col[3].values.all(), 'copy_columns() failed'

# catenate_data()
concat_df = catenate_data([_v, _p])
assert concat_df.shape > _v.shape and concat_df.shape > _p.shape, 'catenate_data() failed'

# Test main()
data, stats, unmapped = __main__(
    network=network,
    transcriptomics_url=transcriptomics_url,
    proteomics_url=proteomics_url,
    metabolomics_url=metabolomics_url)
assert data.shape > _v.shape, '__main__() from prepare_data.py failed'
try:
    data.loc['UBE2N']
except KeyError:
    pass
else:
    raise Exception('__main__() from prepare_data.py failed')
try:
    data.loc['P61088']
except KeyError:
    raise Exception('__main__() from prepare_data.py failed')
try:
    data.loc['ENST00000318066']
except KeyError:
    raise Exception('__main__() from prepare_data.py failed') 
