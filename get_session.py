"""

"""

import sys, os

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

import sqlalchemy.exc
import ConfigParser


def get_session(env):
    """
        looks up the connection parameters in a config file, which must be
        located in the same directory as the python scripts
    """

    filename = 'motifatlas.cfg'
    environments = ['production', 'dev', 'test']
    if env not in environments:
        print 'Error: environment not specified.'
        sys.exit(2)

    script_path = os.path.dirname(os.path.abspath( __file__ ))
    configfile = os.path.join(script_path, filename)
    config = ConfigParser.RawConfigParser()
    config.read(configfile)

    if env == 'production':
        db = 'database'
    elif env == 'test':
        db = 'database_test'
    elif env == 'dev':
        db = 'database_dev'

    engine = create_engine('mysql://' + config.get('database','user')     + ':' +
                                      config.get('database','password') + '@' +
                                      config.get('database','host')     + '/' +
                                      config.get('database', db))


    Session = sessionmaker(bind=engine)
    return Session(), engine
