#!/usr/bin/env python

"""The main script to run BGSU's update pipeline. This takes the name of a
stage to run, options, and then optional list of pdbs. In order for this to run
the configuration file in `conf/motifatlas.json` must exist and be filled out.
There is an example file in `conf/motifatlas.josn.txt`. Notably, a database uri
must be configured. Locations
"""

import os
import sys
import logging
import argparse
from datetime import datetime as dt

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

here = os.path.dirname(__file__)
pymotifs = os.path.abspath(os.path.join(here, '..'))
sys.path.append(pymotifs)

from pymotifs import config
from pymotifs import models
from pymotifs.dispatcher import Dispatcher
from pymotifs.utils import pdb
from pymotifs.utils import known
from pymotifs import email

logger = logging.getLogger(__name__)


def setup_logging(opts):
    log_args = {
        'level': getattr(logging, opts.pop('log_level').upper()),
        'filemode': opts.pop('log_mode'),
        'format': '%(levelname)s:%(asctime)s:%(name)s:%(message)s',
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


def run(name, conf, pdbs, opts):
    setup_logging(opts)

    if opts['before'] or opts['after']:
        pdbs = pdb.RnaPdbsHelper()(dates=(opts['after'], opts['before']))

    if opts.pop('known'):
        pdbs = list(known(config, pdb=False))

    if opts['all']:
        pdbs = pdb.RnaPdbsHelper()()

    if not pdbs:
        logger.warning("Running with no pdb files")

    engine = create_engine(conf['db']['uri'])
    models.reflect(engine)
    Session = sessionmaker(bind=engine)

    dispatcher = Dispatcher(name, conf, Session,
                            skip_dependencies=opts.get('skip_dependencies'),
                            exclude=opts['exclude'])
    dispatcher(pdbs, **opts)


if __name__ == '__main__':
    date = lambda r: dt.strptime(r, '%Y-%m-%d').date()

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('name', metavar='N', help='Name of stage to run')
    parser.add_argument('pdbs', metavar='P', nargs='*', help="PDBs to use")

    parser.add_argument('--config', dest='config',
                        default='conf/motifatlas.json',
                        help="Config file to use")

    # Options about which PDB files to use
    parser.add_argument('--all', dest='all', default=False,
                        action='store_true',
                        help="Use all RNA containing PDBS")
    parser.add_argument('--after-date', dest='after', default=None,
                        type=date, help='Get Pdbs from after the date')
    parser.add_argument('--before-date', dest='before', default=None,
                        type=date, help='Get PDBs from before the date')
    parser.add_argument('--known', action='store_true',
                        help="Use only downloaded pdbs")

    # Options about the type of run
    parser.add_argument('--recalculate', action='store_true',
                        help="Force all data to be recalculated")
    parser.add_argument('--dry-run', action='store_true',
                        help="Do a dry run where we store nothing")
    parser.add_argument('--ignore-time', action='store_true',
                        help='Do not use time for rerunning')

    # Skip options
    parser.add_argument('--skip-dependencies', action='store_true',
                        help='Skip running any dependencies')
    parser.add_argument('--skip-stage', action='append', dest='exclude',
                        help='Name of stage to skip')
    parser.add_argument('--exclude', action='append', dest='skip_pdbs',
                        help='PDB ids to skip')

    parser.add_argument('--no-email', action='store_false',
                        help='Do not send an email')

    # Logging options
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
    opts['exclude'] = set(opts['exclude'] or [])

    conf = config.load(args.config)
    conf['email']['send'] = args.no_email
    run(args.name, conf, args.pdbs, opts)
