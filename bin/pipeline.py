#!/usr/bin/env python

"""The main script to run BGSU's update pipeline. This takes the name of a
stage to run, options, and then optional list of pdbs. In order for this to run
the configuration file in `conf/motifatlas.json` must exist and be filled out.
There is an example file in `conf/motifatlas.josn.txt`. Notably, a database uri
must be configured. Locations
"""

import os
import sys
import argparse
from datetime import datetime as dt

here = os.path.dirname(__file__)
pymotifs = os.path.abspath(os.path.join(here, '..'))
sys.path.append(pymotifs)

from pymotifs import setup


def run(runnable):
    runnable.options.pop('target')
    runnable.dispatcher(runnable.pdbs, **runnable.options)


def inspect(runnable):
    dispatcher = runnable.dispatcher
    ids = runnable.pdbs
    for stage in dispatcher.stages(ids, build=True):
        print(stage.name)


def stages(*args):
    pass


def parser():
    date = lambda r: dt.strptime(r, '%Y-%m-%d').date()
    parser = argparse.ArgumentParser(description=__doc__)
    sps = parser.add_subparsers(title='subcommands',
                                description='known subcommands',
                                help='select a command to run')

    # Add all common options
    parser.add_argument('--config', dest='config',
                        default='conf/motifatlas.json',
                        help="Configuration file to use")
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
    runner = sps.add_parser('run', help='Run a stage')
    runner.set_defaults(target=run)

    inspecter = sps.add_parser('inspect', help='Inspect run')
    inspecter.set_defaults(target=inspect)

    for instance in (runner, inspecter):
        instance.add_argument('name', metavar='N',
                              help='Name of stage to run')

        instance.add_argument('--skip-dependencies', action='store_true',
                              help='Skip running any dependencies')
        instance.add_argument('--skip-stage', action='append',
                              dest='exclude', default=[],
                              help='Name of the stage(s) to skip')
        instance.add_argument('--recalculate', action='append', default=[],
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

    # Setup the bootstrapping parser
    bootstraper = sps.add_parser('bootstrap', help='Fill up a database')
    bootstraper.set_defaults(target=run, stage="update", seed=1,
                             exclude=["units.distances"],
                             config='conf/bootstrap.json')

    # Setup the listing parser
    lister = sps.add_parser('stages', help='List known stages')
    lister.add_argument('pattern', metavar='P', nargs='?',
                        help='Name pattern to use')
    lister.set_defaults(target=stages)

    # Setup the secondary structure import parser
    ss_importer = sps.add_parser('ss.import', help='Import a 2d structure')
    ss_importer.add_argument('filename', metavar='F',
                             help='Secondary structure file to import')
    ss_importer.set_defaults(target=run, name='ss.positions', ids=[])

    # Setup the secondary structure alignment parser
    ss_aligner = sps.add_parser('ss.align', help='Align a PDB, chain to a 2d')
    ss_aligner.add_argument('pdb', metavar='P', help='Name of pdb to use')
    ss_aligner.add_argument('chain', metavar='C', help='Name of the chain')
    ss_aligner.add_argument('file', metavar='S', help='2D structure filename')
    ss_aligner.set_defaults(target=run, name='ss.position_mapping', ids=[])

    # Add the common options
    for instance in [ss_importer, ss_aligner]:
        instance.add_argument('--ss-name', default=None, type=str,
                              help='Name of the 2d diagram')
        instance.add_argument('--recalculate', action='append', default=[],
                              help="Recalculate data for the given stage(s)")

    return parser


def main():
    cli = parser()
    args = cli.parse_args()
    runnable = setup.setup(vars(args))
    args.target(runnable)


if __name__ == '__main__':
    main()
