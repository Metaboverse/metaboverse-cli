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
import argparse
import sklearn
import networkx
import scipy
import numpy
import pandas
import requests
import pickle
import sys
import os

"""Import internal dependencies
"""
try:
    from __init__ import __version__
    from arguments import parse_arguments
    from curate.__main__ import __main__ as curate
    from analyze.__main__ import __main__ as analyze
    from mapper.__main__ import __main__ as mapper
    from target.__main__ import __main__ as curate_target
    from utils import progress_feed, update_session, safestr
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "__version__", os.path.abspath("./metaboverse_cli/__init__.py"))
    version = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(version)
    __version__ = version.__version__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/arguments.py"))
    arguments = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(arguments)
    parse_arguments = arguments.parse_arguments

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/curate/__main__.py"))
    curate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(curate)
    curate = curate.__main__

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/analyze/__main__.py"))
    analyze = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyze)
    analyze = analyze.__main__

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/mapper/__main__.py"))
    mapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mapper)
    mapper = mapper.__main__

    spec = importlib.util.spec_from_file_location(
        "__main__", os.path.abspath("./metaboverse_cli/target/__main__.py"))
    target = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(target)
    curate_target = target.__main__

    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    progress_feed = utils.progress_feed
    update_session = utils.update_session
    safestr = utils.safestr


def main(
        args=None):
    """Run metaboverse-cli
    """

    # Read in arguments
    args, args_dict = parse_arguments(
        args,
        __version__)

    print('Back-end parsed arguments:')
    print(args_dict)
    print()

    if args_dict['cmd'] == 'metaboliteMapper':

        print('Generating metabolite mapper...')
        mapper(args_dict)
        sys.exit(1)

    # Run metaboverse-curate
    elif args_dict['cmd'] == 'curate' or args_dict['cmd'] == 'electrum':

        if args_dict['cmd'] == 'curate':
            print('Generating Metaboverse-compatible database...')
        elif args_dict['cmd'] == 'electrum':
            print('Generating Electrum-compatible database...')

        if safestr(args_dict['organism_curation']) != 'None':
            progress_feed(
                args_dict=args_dict,
                process="curate",
                amount=50)
            # Update args_dict with path for network model
            with open(args_dict['organism_curation'], 'rb') as network_file:
                network = pickle.load(network_file)
                args_dict['organism_id'] = network['organism_id']
                if args_dict['output_file'] == None \
                        or args_dict['output_file'] == "None" \
                        or args_dict['output_file'] == "find":
                    args_dict['output_file'] = args_dict['output'] \
                        + args_dict['organism_id'] \
                        + '.mvrs'
                args_dict['network'] = args_dict['organism_curation']

                # add variables back to session data json file
                session_file = args_dict['session_data']
                update_session(
                    session_file=session_file,
                    key='organism_id',
                    value=args_dict['organism_id'])
                update_session(
                    session_file=session_file,
                    key='output_file',
                    value=args_dict['output_file'])
                update_session(
                    session_file=session_file,
                    key='database_url',
                    value=args_dict['output_file'])
            print('Skipping organism network modeling as one was provided by'
                  + ' the user...')

        else:
            print('Curating network model...')
            if 'model_file' in args_dict \
                    and safestr(args_dict['model_file']) == 'None':
                args_dict['model_file'] = args_dict['output'] \
                    + args_dict['organism_id'] \
                    + '.mvdb'

            args_dict['network'] = args_dict['model_file']
            args_dict = curate(args_dict)
            # sys.stdout.flush()

        if 'output_file' in args_dict \
                and safestr(args_dict['output_file']) == 'None':
            args_dict['output_file'] = args_dict['output'] \
                + args_dict['organism_id'] \
                + '.mvrs'

        print('Curating data onto the network model...')
        if args_dict['cmd'] == 'curate':
            if 'output_file' in args_dict \
                    and safestr(args_dict['output_file']) == 'None':
                args_dict['output_file'] = args_dict['output'] \
                    + args_dict['organism_id'] \
                    + '.mvrs'
            analyze(args_dict)
        elif args_dict['cmd'] == 'electrum':
            if 'output_file' in args_dict \
                    and safestr(args_dict['output_file']) == 'None':
                args_dict['output_file'] = args_dict['output'] \
                    + args_dict['organism_id'] \
                    + '-latest.eldb'
            curate_target(args_dict)

        sys.exit(1)

    # Print some error messaging
    else:
        raise Exception('Invalid sub-module selected')


"""Run main
"""
if __name__ == '__main__':

    sys.exit(main() or 0)
