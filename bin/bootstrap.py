#!/usr/bin/env python

"""
This is a script to bootstrap up a testing database. This script runs stages of
the pipeline on a specific subset of structures. This reads the config file in
`conf/bootstrap.json` for the database and pdb files to bootstrap. For an
example file see `conf/bootstrap.json.txt`. Note that the bootstrapping process
is not perfect, if matlab is not available then several stages, such as
importing interactions, will be skipped. Also, this requires that a database
already be created with the correct schema. This can be achieved by using the
migration scripts in our database-migrations repo.
"""

import argparse

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import pipeline as pipe

from pymotifs import models
from pymotifs.dispatcher import Dispatcher


def main(config, opts):
    pipe.setup_logging(opts)

    engine = create_engine(config['db']['uri'])
    models.reflect(engine)
    Session = sessionmaker(bind=engine)

    dispatcher = Dispatcher("update", config, Session)
    dispatcher(config['pdbs'], **opts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', dest='config',
                        default='conf/bootstrap.json',
                        help="Config file to use")
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
        if arg != 'config':
            opts[arg] = value

    config = pipe.load_config(args.config)
    main(config, opts)
