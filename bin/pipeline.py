#!/usr/bin/env python

"""The main script to run BGSU's update pipeline. This takes the name of a
stage to run, options, and then optional list of pdbs. In order for this to run
the configuration file in `conf/motifatlas.json` must exist and be filled out.
There is an example file in `conf/motifatlas.josn.txt`. Notably, a database uri
must be configured. Locations
"""

import os
import sys
import random
import logging
import argparse
from datetime import datetime as dt
from copy import deepcopy

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

here = os.path.dirname(__file__)
pymotifs = os.path.abspath(os.path.join(here, '..'))
sys.path.append(pymotifs)

import pymotifs as pym


def setup_logging(options):
    opts = deepcopy(options)
    log_args = {
        'level': getattr(logging, opts.pop('log_level').upper()),
        'filemode': opts.pop('log_mode'),
        'format': '%(levelname)s:%(name)s:%(asctime)s:%(message)s',
    }

    filename = opts.pop('log_file')
    if filename:
        log_args['filename'] = filename

    logging.basicConfig(**log_args)
    # logging.captureWarnings(True)
    base = logging.getLogger()
    pool_logger = logging.getLogger('sqlalchemy.pool')
    pool_logger.setLevel(logging.ERROR)
    for handler in base.handlers:
        pool_logger.addHandler(handler)

    return opts


def setup(args):
    options = {
        'exclude': set()
    }

    stage = args.name
    ids = args.ids
    config = pym.config.load(args.config)
    arguments = set(['ids', 'stage', 'config'])
    for arg, value in vars(args).items():
        if arg not in arguments:
            options[arg] = value

    opts = setup_logging(options)
    engine = create_engine(config['db']['uri'])
    pym.models.reflect(engine)
    Session = sessionmaker(bind=engine)
    opts['dispatcher'] = pym.dispatcher.Dispatcher(stage, config, Session,
                                                   **options)

    if not ids:
        if options['before'] or options['after']:
            dates = (options['after'], options['before'])
            ids = pym.utils.pdb.RnaPdbsHelper()(dates=dates)

        if options.pop('known'):
            ids = list(pym.utils.known(config, pdb=False))

        if options['all']:
            ids = pym.utils.pdb.RnaPdbsHelper()()

    if 'seed' in options:
        random.seed(options.pop('seed'))

    return ids, opts


def run(ids, options):
    dispatcher = options.pop('dispatcher')
    dispatcher(ids, **options)


def inspect(ids, options):
    dispatcher = options.pop('dispatcher')
    for stage in dispatcher.stages(ids, build=True):
        print(stage.name)


def stages(*args):
    pass


def parser():
    date = lambda r: dt.strptime(r, '%Y-%m-%d').date()
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='known subcommands',
                                       help='select a command to run')

    # Add all common options
    parser.add_argument('--log-file', dest='log_file', default='',
                        help="Log file to use")
    parser.add_argument('--log-level', dest='log_level',
                        default='debug',
                        choices=['debug', 'info', 'warning', 'error'],
                        help="Logging level to use")
    parser.add_argument('--log-mode', dest='log_mode', default='a',
                        choices=['w', 'a'],
                        help='Mode to open the  logging file')

    # Setup the common parts of the run and inspect commands
    runner = subparsers.add_parser('run', help='Run a stage')
    runner.set_defaults(target=run)

    inspecter = subparsers.add_parser('inspect', help='Inspect run')
    inspecter.set_defaults(target=inspect)

    for instance in (runner, inspecter):
        instance.add_argument('name', metavar='N', nargs=1,
                              help='Name of stage to run')

        instance.add_argument('--skip-dependencies', action='store_true',
                              help='Skip running any dependencies')
        instance.add_argument('--skip-stage', action='append',
                              dest='exclude',
                              help='Name of the stage(s) to skip')
        instance.add_argument('--recalculate', action='store_true',
                              help="Recalculate data for the given stage(s)")
        instance.add_argument('--all', dest='all', default=False,
                              action='store_true',
                              help="Use all RNA containing PDBS")
        instance.add_argument('--known', action='store_true',
                              help="Use only downloaded cif files")
        instance.add_argument('--after-date', dest='after', default=None,
                              type=date, help='Get Pdbs from after the date')
        instance.add_argument('--before-date', dest='before', default=None,
                              type=date, help='Get PDBs from before the date')
        instance.add_argument('--exclude', action='append', dest='skip_pdbs',
                              help='PDB id(s) to skip')
        instance.add_argument('--ignore-time', action='store_true',
                              help='Do not use time for rerunning')

    # Setup the run parser
    runner.add_argument('ids', metavar='P', nargs='*', help="PDBs to use")
    runner.add_argument('--dry-run', action='store_true',
                        help="Do a dry run where we alter nothing")

    bootstraper = subparsers.add_parser('bootstrap', help='fill a database')
    bootstraper.set_defaults(target=run, stage="update", seed=1,
                             exclude=["units.distances"],
                             config='conf/bootstrap.json')

    lister = subparsers.add_parser('stages', help='List known stages')
    lister.add_argument('pattern', metavar='P', nargs='?',
                        help='Name pattern to use')
    lister.set_defaults(target=stages)

    return parser


def main():
    cli = parser()
    args = cli.parse_args()
    ids, options = setup(args)
    args.target(ids, options)


if __name__ == '__main__':
    main()
