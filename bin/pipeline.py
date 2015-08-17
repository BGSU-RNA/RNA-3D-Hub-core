#!/usr/bin/env python

import os
import sys
import json
import logging
import argparse
import collections

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

here = os.path.dirname(__file__)
pymotifs = os.path.abspath(os.path.join(here, '..'))
sys.path.append(pymotifs)

from pymotifs import models
from pymotifs.dispatcher import Dispatcher
from pymotifs.utils import pdb
from pymotifs.utils import known
from pymotifs import email

logger = logging.getLogger(__name__)

DESCRIPTION = """
The main script to run BGSU's update pipeline. This takes the name of a stage
to run, a optional list of pdb and then options.
"""


def setup_logging(opts):
    log_args = {
        'level': getattr(logging, opts.pop('log_level').upper()),
        'filemode': opts.pop('log_mode')
    }
    if opts['log_file']:
        log_args['filename'] = opts.pop('log_file')

    logging.basicConfig(**log_args)
    # logging.captureWarnings(True)
    base = logging.getLogger()
    pool_logger = logging.getLogger('sqlalchemy.pool')
    pool_logger.setLevel(logging.ERROR)
    for handler in base.handlers:
        pool_logger.addHandler(handler)


def load_config(filename):
    config = collections.defaultdict(dict)
    with open(filename, 'rb') as raw:
        config.update(json.load(raw))
        return config


def run(name, config, pdbs, opts):
    setup_logging(opts)

    if opts['all']:
        pdbs = pdb.RnaPdbsHelper()()

    if opts.pop('known'):
        pdbs = list(known(config, pdb=False))

    if not pdbs:
        logger.warning("Running with no pdb files")

    engine = create_engine(config['db']['uri'])
    models.reflect(engine)
    Session = sessionmaker(bind=engine)

    dispatcher = Dispatcher(name, config, Session,
                            skip_dependencies=opts.get('skip_dependencies'))

    try:
        dispatcher(pdbs, **opts)
        if opts.get('no_email', True):
            email.send(name, log_file=opts.get('log_file'))
    except Exception as err:
        if opts.get('no_email', True):
            email.send(name, log_file=opts.get('log_file'), error=err)
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument('name', metavar='N', help='Name of stage to run')
    parser.add_argument('pdbs', metavar='P', nargs='*', help="PDBs to use")

    parser.add_argument('--all', dest='all', default=False,
                        action='store_true',
                        help="Use all RNA containing PDBS")
    parser.add_argument('--config', dest='config',
                        default='conf/motifatlas.json',
                        help="Config file to use")
    parser.add_argument('--recalculate', action='store_true',
                        help="Force all data to be recalculated")
    parser.add_argument('--known', action='store_true',
                        help="Use only downloaded pdbs")
    parser.add_argument('--dry-run', action='store_true',
                        help="Do a dry run where we store nothing")
    parser.add_argument('--skip-dependencies', action='store_true',
                        help='Skip running any dependencies')

    parser.add_argument('--no-email', action='store_true',
                        help='Do not send an email')

    parser.add_argument('--log-file', dest='log_file', default='',
                        help="Log file to use")
    parser.add_argument('--log-level', dest='log_level', default='debug',
                        choices=['debug', 'info', 'warning', 'error'],
                        help="Logging level to use")
    parser.add_argument('--log-mode', dest='log_mode', default='a',
                        choices=['w', 'a'],
                        help='Mode to open the  logging file')

    args = parser.parse_args()
    opts = {}
    for arg, value in vars(args).items():
        if arg != 'pdbs' and arg != 'name' and arg != 'config':
            opts[arg] = value

    config = load_config(args.config)
    run(args.name, config, args.pdbs, opts)
