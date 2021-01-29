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
import re
import sys
import os

def read_data(
        url,
        delimiter='\t'):
    """Expected to contain an even number of columns for the data type
    """

    data = pd.read_csv(
        url,
        sep=delimiter,
        index_col=0)

    if len(data.columns.tolist()) % 2 != 0:
        raise Exception('Improperly formatted datatable provided: ', url)

    return data


def format_data(
        data,
        reference):
    """Format data for processing
    0) Transcriptomics
    - Accepts entities mapping to ensembl genes and their synonyms
    - Expects fold change values and their p-values side-by-side
    - Columns should be ordered in intended order (for multi-conditions)
    - gene_dictionary is an Ensembl GTF file
    - Columns should be condition_fc, condition_p, etc.
    1) Proteomics
    - Accepts entities mapping to uniprot proteins and their synonyms
    - Expects fold change values and their p-values side-by-side
    - Columns should be ordered in intended order (for multi-conditions)
    - gene_dictionary is an Ensembl GTF file
    - Columns should be condition_fc, condition_p, etc.
    """

    data_output = data.copy()
    data_output.index = data_output.index.str.upper()

    reference_ids = list(reference.values())

    data_output.index = data_output.index.to_series().replace(reference)
    data_unmapped = data_output.copy()

    data_unmapped = data_unmapped[~data_unmapped.index.isin(reference_ids)]

    return data_output, data_unmapped


def output_unmapped(
        data,
        url,
        delimiter='\t'):
    """Output unmapped entities for user information
    """

    if len(data.index.tolist()):
        data.to_csv(
            url[:-4] + '_unmapped.txt',
            sep=delimiter)


def extract_data(
        data):
    """Seperate out fold change and p-values
    - Formatting follows fold change, p-values, etc...
    """
    data_c = data.copy()
    data_c = data_c.dropna(axis=0)
    data_c.index = [d.lstrip().rstrip() for d in data_c.index.tolist()]

    _values = data_c.T[::2].T
    _stats = data_c.T[1::2].T

    _values.columns = [x for x in range(len(_values.columns))]
    _stats.columns = [x for x in range(len(_stats.columns))]

    return _values, _stats


def broadcast_transcriptomics(
        transcriptomics,
        transcriptomics_stats,
        gene_dictionary,
        protein_dictionary):
    """If transcript levels are provided and no protein levels are given,
    the user can opt to broadcast the gene expression levels to the protein
    nodes. If an enzyme node is constituted of several entities, the user can
    choose to use the average or use the most extreme value (the fc furthest
    from 0)
    """

    proteomics = transcriptomics.copy()
    proteomics_stats = transcriptomics_stats.copy()

    uniprot_dict = {}
    for x in protein_dictionary.keys():
        uniprot_dict[protein_dictionary[x]] = x

    proteomics.index = proteomics.index.to_series().replace(
        gene_dictionary)
    proteomics_stats.index = proteomics_stats.index.to_series().replace(
        gene_dictionary)

    proteomics.index = proteomics.index.to_series().replace(
        uniprot_dict)
    proteomics_stats.index = proteomics_stats.index.to_series().replace(
        uniprot_dict)

    return proteomics, proteomics_stats


def copy_columns(
        data,
        stats,
        _max):
    """Copy data columns for every time point provided by the user
    """

    data_c = data.copy()
    stats_c = stats.copy()

    counter = 1

    while counter < _max:
        data_c[counter] = data_c[0]
        stats_c[counter] = stats_c[0]
        counter += 1

    return data_c, stats_c


def catenate_data(
        array):
    """Combine data
    - average duplicates
    - remove NAs
    Return: combined dataframe where indices are species IDs
    """

    combined = pd.concat(array)
    combined = combined.dropna(axis=0)
    combined = combined.sort_index()

    removers = []  # Remove non-numbers
    for idx, row in combined.iterrows():
        for x in row:
            try:
                float(x)
            except:
                removers.append(idx)

    combined = combined[~combined.index.isin(removers)]

    for col in combined.columns.tolist():
        combined[col] = combined[col].astype(float)

    combined = combined.groupby(combined.index).mean()

    return combined


def __main__(
        network,
        transcriptomics_url,
        proteomics_url,
        metabolomics_url,
        database_source='reactome'):
    """Get user data and preprocess
    """

    # Process transcriptomics
    if transcriptomics_url.lower() != 'none':

        transcriptomics = read_data(
            url=transcriptomics_url)
        e_sym = {}
        if database_source.lower() == 'reactome':
            for k, v in network['ensembl_synonyms'].items():
                if 'phospho-' in v and '-phospho-' not in v:
                    v = v.replace('phospho-', '')
                if '(' in v and ')' in v:
                    v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
                if '  ' in v:
                    v = v.replace('  ', ' ')
                v = v.strip()
                e_sym[v.upper()] = k
                e_sym[k.upper()] = k
            transcriptomics, transcriptomics_unmapped = format_data(
                data=transcriptomics,
                reference=e_sym)
            output_unmapped(
                data=transcriptomics_unmapped,
                url=transcriptomics_url)
        transcriptomics, transcriptomics_stats = extract_data(
            data=transcriptomics)
        transcriptomics_length = len(transcriptomics.columns.tolist())
    else:
        transcriptomics_length = 0

    # Process proteomics
    if proteomics_url.lower() != 'none':

        proteomics = read_data(
            url=proteomics_url)
        u_sym = {}
        if database_source.lower() == 'reactome':
            for k, v in network['uniprot_synonyms'].items():
                if 'phospho-' in v and '-phospho-' not in v:
                    v = v.replace('phospho-', '')
                if '(' in v and ')' in v:
                    v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
                if '  ' in v:
                    v = v.replace('  ', ' ')
                v = v.strip()
                u_sym[v.upper()] = k
                u_sym[k.upper()] = k
            proteomics, proteomics_unmapped = format_data(
                data=proteomics,
                reference=u_sym)
            output_unmapped(
                data=proteomics_unmapped,
                url=proteomics_url)
        proteomics, proteomics_stats = extract_data(
            data=proteomics)
        proteomics_length = len(proteomics.columns.tolist())
    else:
        proteomics_length = 0

    # Process metabolomics
    if metabolomics_url.lower() != 'none':

        metabolomics = read_data(
            url=metabolomics_url)
        metabolomics, metabolomics_stats = extract_data(
            data=metabolomics)
        metabolomics_length = len(metabolomics.columns.tolist())
    else:
        metabolomics_length = 0

    # Check for broadcasting
    if proteomics_url.lower() == 'none' \
            and transcriptomics_url.lower() != 'none':

        proteomics, proteomics_stats = broadcast_transcriptomics(
            transcriptomics=transcriptomics,
            transcriptomics_stats=transcriptomics_stats,
            gene_dictionary=network['ensembl_synonyms'],
            protein_dictionary=network['uniprot_synonyms'])

    # Allow for unequal filling
    lengths = []
    for x in [
        transcriptomics_length,
        proteomics_length,
        metabolomics_length,
    ]:
        if x != 0:
            lengths.append(x)
        else:
            lengths.append(1)

    if sum(lengths) != 3 and len(list(set(lengths))) == 2:

        _max = max(lengths)

        if transcriptomics_url.lower() != 'none' \
                and len(transcriptomics.columns.tolist()) != _max \
                and len(transcriptomics.columns.tolist()) == 1:
            transcriptomics, transcriptomics_stats = copy_columns(
                data=transcriptomics,
                stats=transcriptomics_stats,
                _max=_max)

        if proteomics_url.lower() != 'none' \
                and len(proteomics.columns.tolist()) != _max \
                and len(proteomics.columns.tolist()) == 1:
            proteomics, proteomics_stats = copy_columns(
                data=proteomics,
                stats=proteomics_stats,
                _max=_max)

        if metabolomics_url.lower() != 'none' \
                and len(metabolomics.columns.tolist()) != _max \
                and len(metabolomics.columns.tolist()) == 1:
            metabolomics, metabolomics_stats = copy_columns(
                data=metabolomics,
                stats=metabolomics_stats,
                _max=_max)

    elif len(list(set(lengths))) == 1:
        pass
    else:
        print("When providing multi-omic timecourse data with unequals times, other omics types must match the maximum number of time points or must only provide one time point.")

    # Initialize array of filled data
    data_array = []
    stats_array = []

    if transcriptomics_url.lower() != 'none':
        data_array.append(transcriptomics)
        stats_array.append(transcriptomics_stats)
    if proteomics_url.lower() != 'none':
        data_array.append(proteomics)
        stats_array.append(proteomics_stats)
    if metabolomics_url.lower() != 'none':
        data_array.append(metabolomics)
        stats_array.append(metabolomics_stats)

    # Combine data
    data = catenate_data(
        array=data_array)

    stats = catenate_data(
        array=stats_array)

    unmapped = {}
    try:
        unmapped['transcriptomics_unmapped'] = \
            transcriptomics_unmapped.index.tolist()
    except:
        unmapped['transcriptomics_unmapped'] = []
    try:
        unmapped['proteomics_unmapped'] = \
            proteomics_unmapped.index.tolist()
    except:
        unmapped['proteomics_unmapped'] = []

    return data, stats, unmapped
