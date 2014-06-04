import re
import os
import ConfigParser
import logging

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base


def camelize_classname(base, tablename, table):
    """Produce a 'camelized' class name, e.g.
    'words_and_underscores' -> 'WordsAndUnderscores'
    """

    return str(tablename[0].upper() +
               re.sub(r'_(\w)', lambda m: m.group(1).upper(), tablename[1:]))


def get_engine(filename='motifatlas.cfg'):
    """
        looks up the connection parameters in a config file, which must be
        located in the same directory as the python scripts.
    """

    script_path = os.path.dirname(os.path.abspath(__file__))
    configfile = os.path.join(script_path, filename)
    config = ConfigParser.RawConfigParser()
    config.read(configfile)

    uri = config.get('database', 'uri')
    logging.info('Connecting to the `%s` database', uri)

    return create_engine(uri)


def tables(engine):
    Base = automap_base()
    Base.prepare(engine, reflect=True, classname_for_table=camelize_classname)
    for klass in Base.classes:
        yield klass.__name__, klass


engine = get_engine()
Session = sessionmaker(bind=engine)
session = Session()

for name, klass in tables(engine):
    globals()[name] = klass
