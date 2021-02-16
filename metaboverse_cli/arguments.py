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
from textwrap import dedent
import argparse
import sys
import os

"""Import internal dependencies
"""
try:
    from __init__ import __version__
    from utils import check_directories
    from utils import check_curate
    from utils import argument_checks
    from utils import safestr
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "__version__", os.path.abspath("./metaboverse_cli/__init__.py"))
    version = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(version)
    __version__ = version.__version__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    check_directories = utils.check_directories
    check_curate = utils.check_curate
    argument_checks = utils.argument_checks
    safestr = utils.safestr

__path__ = os.path.dirname(os.path.realpath(__file__))
url = 'https://raw.githubusercontent.com/j-berg/Metaboverse/master/metaboverse/__init__.py'

description_table = """\
    The Metaboverse sub-modules can be accessed by executing:
        'metaboverse sub-module_name arg1 arg2 ...'
    Sub-module help can be displayed by executing:
    'metaboverse sub-module_name --help'
    Sub-module descriptions:
        +-----------------------+--------------------------------------------+
        |   metaboliteMapper    |   Process a metabolite mapper database     |
        +-----------------------+--------------------------------------------+
        |   curate              |   Curate network with optional user data   |
        +-----------------------+--------------------------------------------+
"""

def check_arguments(
        args_dict):
    """Check arguments
    """

    # Run general checks
    args_dict = argument_checks(args_dict)

    # Run sub-module specific checks
    if args_dict['cmd'] == 'curate':
        check_curate(args_dict)
    elif args_dict['cmd'] == 'metaboliteMapper':
        pass
    elif args_dict['cmd'] == 'electrum':
        pass
    else:
        raise Exception('Invalid sub-module selected')

    if 'output' not in args_dict \
            or args_dict['output'] == None \
            or args_dict['output'].lower() == 'none':
        args_dict['output'] = args_dict['output_file'].rsplit(os.path.sep, 1)[
            0]

    if not args_dict['output'].endswith(os.path.sep):
        args_dict['output'] = args_dict['output'] + os.path.sep

    args_dict['output'] = check_directories(
        input=args_dict['output'],
        argument='output')

    # Print out user commands to log file
    print('Metaboverse-CLI version: ' + safestr(__version__))
    print('======================\nUser commands summary:\n======================')

    for key, value in args_dict.items():
        print(safestr(key) + ': ' + safestr(value))
    print('=====================\nEnd commands summary\n=====================\n')

    return args_dict


"""Parse arguments
Will print help menu if no arguments are provided
"""


def parse_arguments(
        args,
        __version__):

    # Require user input
    if args is None:
        args = sys.argv[1:]

    # Initialize main parser
    parser = argparse.ArgumentParser(
        prog='metaboverse',
        description=dedent(description_table),
        formatter_class=argparse.RawTextHelpFormatter)

    # Optional main arguments
    parser.add_argument(
        '-v', '--version',
        help='Print installed version to stout',
        action='version',
        version='%(prog)s ' + safestr(__version__))

    # Add sub-parsers
    subparser = parser.add_subparsers(dest='cmd')

    # electrum parser
    electrum_parser = subparser.add_parser(
        'electrum',
        description='Curate Electrum-compatible database',
        add_help=False)

    # electrum required arguments
    electrum_reqs = electrum_parser.add_argument_group('required arguments')
    electrum_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory (default: current working directory)',
        metavar='<path>',
        type=str,
        required=True)
    electrum_reqs.add_argument(
        '-d', '--data',
        help='Path and filename of MIDAS interaction database',
        metavar='<path/filename.txt>',
        type=str,
        required=True)
    electrum_reqs.add_argument(
        '-s', '--organism_id',
        help='Reactome species ID',
        metavar='<organism_id>',
        type=str,
        default='HSA',
        required=True)

    # electrum optional arguments
    electrum_opts = electrum_parser.add_argument_group('optional arguments')
    electrum_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    electrum_opts.add_argument(
        '--database_source',
        help='A string indicating the database source (default: "reactome"; other options: "biomodels/bigg")',
        metavar='<source_name>',
        type=str,
        default='reactome',
        required=False)
    electrum_opts.add_argument(
        '-c', '--organism_curation',
        help='Path and name for organism curation file',
        metavar='<path/filename.mvdb; path/filename.xml; path/filename.sbml>',
        type=str,
        required=False)
    electrum_opts.add_argument(
        '-i', '--model_file',
        help='Path and name for organism curation file output. If --organism_curation is used, this argument will be ignored.',
        metavar='<path/filename.mvdb>',
        type=str,
        required=False)
    electrum_opts.add_argument(
        '-f', '--output_file',
        help='Path and name for output database file',
        metavar='<path/filename.mvrs>',
        type=str,
        required=False)

    # metaboliteMapper parser
    mapper_parser = subparser.add_parser(
        'metaboliteMapper',
        description='Curate metabolite mapper',
        add_help=False)

    # metaboliteMapper required arguments
    mapper_reqs = mapper_parser.add_argument_group('required arguments')
    mapper_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory (default: current working directory)',
        metavar='<path>',
        type=str,
        required=True)

    # Curate parser
    curate_parser = subparser.add_parser(
        'curate',
        description='Curate biological network',
        add_help=False)

    # Curate required arguments
    curate_reqs = curate_parser.add_argument_group('required arguments')
    curate_reqs.add_argument(
        '-o', '--output',
        help='Path to output directory (default: current working directory)',
        metavar='<path>',
        type=str,
        required=True)
    curate_reqs.add_argument(
        '-s', '--organism_id',
        help='Reactome species ID',
        metavar='<organism_id>',
        type=str,
        default='HSA',
        required=True)

    # Curate optional arguments
    curate_opts = curate_parser.add_argument_group('optional arguments')
    curate_opts.add_argument(
        '-h', '--help',
        action='help',
        help='Show help message and exit')
    curate_opts.add_argument(
        '--organism_curation_file',
        help='Path and name for organism curation file output. If --organism_curation is used, this argument will be ignored.',
        metavar='<path/filename.mvdb>',
        type=str,
        required=False)
    curate_opts.add_argument(
        '--neighbor_dictionary_file',
        help='Path and name for organism neighbor dictionary file. If not provided, Metaboverse will search the SourceForge archives for a compatible version.',
        metavar='<path/filename.nbdb>',
        type=str,
        required=False)
    curate_opts.add_argument(
        '--graph_template_file',
        help='Path and name for organism graph template file. If not provided, Metaboverse will search the SourceForge archives for a compatible version.',
        metavar='<path/filename_template.mvrs>',
        type=str,
        required=False)
    curate_opts.add_argument(
        '-f', '--output_file',
        help='Path and name for output database file',
        metavar='<path/filename.mvrs>',
        type=str,
        required=False)
    curate_opts.add_argument(
        '-t', '--transcriptomics',
        help='Path and filename of RNA-Seq data - refer to documentation for details on formatting and normalization',
        metavar='<path/filename>',
        type=str,
        default='None',
        required=False)
    curate_opts.add_argument(
        '-p', '--proteomics',
        help='Path and filename of proteomics data - refer to documentation for details on formatting and normalization',
        metavar='<path/filename>',
        type=str,
        default='None',
        required=False)
    curate_opts.add_argument(
        '-m', '--metabolomics',
        help='Path and filename of metabolomics data - refer to documentation for details on formatting and normalization',
        metavar='<path/filename>',
        type=str,
        default='None',
        required=False)
    curate_opts.add_argument(
        '--database_source',
        help='A string indicating the database source (default: "reactome"; other options: "biomodels/bigg")',
        metavar='<source_name>',
        type=str,
        default='reactome',
        required=False)
    # curate_opts.add_argument(
    #    '-a', '--additional_reactions',
    #    help = 'Path and filename of additional reaction table. See #documentation for more details on appropriate file formatting.',
    #    metavar = '<file.txt, file.tsv>',
    #    type = str,
    #    default = 'None',
    #    required = False)
    curate_opts.add_argument(
        '--force_new_curation',
        help='Force all intermediate database files to be freshly created.',
        action='store_true',
        required=False)
    curate_opts.add_argument(
        '--collapse_with_modifiers',
        help='Include modifiers when considering a potential reaction collapse.',
        action='store_true',
        required=False)
    curate_opts.add_argument(
        '--broadcast_genes',
        help='Broadcast gene values to proteins where protein values not available.',
        action='store_true',
        required=False)
    curate_opts.add_argument(
        '--broadcast_metabolites',
        help='Broadcast metabolites values to protein complexes.',
        action='store_true',
        required=False)
    curate_opts.add_argument(
        '--experiment_type',
        help='Specify experiment type',
        metavar='<default/timecourse/flux/multi-condition>',
        type=str,
        required=False)
    curate_opts.add_argument(
        '--experiment_name',
        help='Specify experiment name',
        metavar='<default/timecourse/flux/multi-condition>',
        type=str,
        required=False)
    curate_opts.add_argument(
        '--labels',
        help='Comma separated list of labels to use for multi-condition or timecourse experiments',
        metavar='cond1, cond2, cond3, ...',
        type=str,
        required=False)
    curate_opts.add_argument(
        '--blocklist',
        help='Comma separated list of names for metabolites, etc., to ignore in the analysis and visualization.',
        metavar='H+, H2O, etc...',
        type=str,
        required=False)
    curate_opts.add_argument(
        '--session_data',
        help='Path and filename to session data file',
        metavar='<path/filename>',
        type=str,
        required=False)
    curate_opts.add_argument(
        '--progress_log',
        help='Path and filename to progress log file',
        metavar='<path/filename>',
        type=str,
        required=False)

    # Get arguments are print help if no arguments provided
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse arguments into NameSpace
    args = parser.parse_args(args)

    # Collect subargs and package, add metaboverse script path to parameter dictionary
    args_dict = vars(args)
    args_dict['path'] = safestr(__path__)
    if not args_dict['path'].endswith(os.path.sep):
        args_dict['path'] = args_dict['path'] + os.path.sep

    # Check inputs validity
    args_dict = check_arguments(
        args_dict)

    print('Back-end parsed arguments:')
    print(args_dict)
    print()

    return args, args_dict
