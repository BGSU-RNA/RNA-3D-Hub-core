#!/usr/bin/env python

"""
A script to generate a report for the nr set. This basically computes what the
nr set would be without writing to the database.
"""

import os
import csv
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
from pymotifs.utils import pdb
from pymotifs.utils import known

logger = logging.getLogger(__name__)


def setup_logging(opts):
    log_args = {
        'level': getattr(logging, opts.pop('log_level').upper()),
        'filemode': opts.pop('log_mode')
    }
    if opts['log_file']:
        log_args['filename'] = opts.pop('log_file')

    logging.basicConfig(**log_args)
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


def format_group(out, groups):
    headers = ['group', 'id', 'length', 'source', 'sequence',
               'discrepancy']

    def as_row(index, data):
        row = dict(data)
        # print(data)
        row['group'] = 'group-%i' % index
        # row['id'] = group['id']
        row['length'] = row['chains'][0]['exp_length']
        row['source'] = ','.join(str(c['source']) for c in row['chains'])
        row['discrepancy'] = ''
        row['sequence'] = row['chains'][0]['sequence']
        return row

    writer = csv.DictWriter(out, headers, delimiter="\t",
                            extrasaction="ignore")

    writer.writerow(dict(zip(headers, headers)))
    for index, group in enumerate(groups):
        writer.writerow(as_row(index, group['representative']))
        for member in group['members']:
            writer.writerow(as_row(index, member))


def format_compare(out, groups):
    grouped = []
    for group in groups:
        chains = [c['id'].split(',')[0] for c in group['members']]
        chains.append(group['representative']['id'].split(',')[0])
        grouped.append(sorted(chains))

    for entry in sorted(grouped):
        out.write(" ".join(entry))
        out.write("\n")


def save(config, filename, groups):
    with open(filename, 'wb') as raw:
        if config['format'] == 'group':
            return format_group(raw, groups)
        elif config['format'] == 'compare':
            return format_compare(raw, groups)

        sys.stderr.write("Bad format\n")
        sys.exit(1)


def run(config, filename, pdbs, opts):
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

    from pymotifs.nr.groups import simplified
    grouper = simplified.Grouper(config, Session)
    groups = grouper(pdbs, **opts)
    save(opts, filename, groups)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('filename', metavar='filename',
                        help="File to write groups to")
    parser.add_argument('pdbs', metavar='P', nargs='*', help="PDBs to use")

    parser.add_argument('--all', dest='all', default=False,
                        action='store_true',
                        help="Use all RNA containing PDBS")
    parser.add_argument('--config', dest='config',
                        default='conf/motifatlas.json',
                        help="Config file to use")
    parser.add_argument('--known', action='store_true',
                        help="Use only downloaded pdbs")
    parser.add_argument('--format', dest='format', default='group',
                        choices=['group', 'compare'])

    parser.add_argument('--log-file', dest='log_file', default='',
                        help="Log file to use")
    parser.add_argument('--log-level', dest='log_level', default='debug',
                        choices=['debug', 'info', 'warning', 'error'],
                        help="Logging level to use")
    parser.add_argument('--log-mode', dest='log_mode', default='w',
                        choices=['w', 'a'],
                        help='Mode to open the  logging file')

    args = parser.parse_args()
    opts = {}
    for arg, value in vars(args).items():
        if arg != 'pdbs' and arg != 'name' and arg != 'config':
            opts[arg] = value

    config = load_config(args.config)
    run(config, args.filename, args.pdbs, opts)
