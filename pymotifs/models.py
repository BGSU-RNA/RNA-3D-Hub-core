import re
import logging

from sqlalchemy import Table
from sqlalchemy import Column
from sqlalchemy import String
from sqlalchemy import Integer
from sqlalchemy import MetaData

from sqlalchemy.ext.declarative import declarative_base


logger = logging.getLogger(__name__)

metadata = MetaData()
Base = declarative_base(metadata=metadata)


class TempPdbs(Base):
    """A class used to when storing pdb_ids temporarily
    """
    __tablename__ = 'temp_pdbs'
    __table_args__ = {'prefixes': ['TEMPORARY']}
    pdb_id = Column(String(4), primary_key=True)


def camelize_classname(tablename):
    return str(tablename[0].upper() +
               re.sub(r'_(\w)', lambda m: m.group(1).upper(), tablename[1:]))


def should_reflect(tablename, *args):
    name = camelize_classname(tablename)
    return name not in globals()


def define_missing_views(metadata):
    """A functionn to define the primary keys for views. We need to do this so
    we can reflect these views.
    """

    Table('exp_seq_pdb', metadata,
          Column('exp_seq_id', String, primary_key=True),
          Column('pdb_id', String, primary_key=True),
          Column('chain_id', Integer, primary_key=True),
          extend_existing=True)

    Table('correspondence_pdbs', metadata,
          Column('correspondence_id', String, primary_key=True),
          Column('chain_id_1', Integer, primary_key=True),
          Column('chain_id_2', Integer, primary_key=True),
          extend_existing=True)

    Table('correspondence_units', metadata,
          Column('correspondence_id', String, primary_key=True),
          Column('unit_id_1', String, primary_key=True),
          Column('unit_id_2', String, primary_key=True),
          extend_existing=True)


def reflect(engine):
    """Reflect all tables/views from the database into python. This cannot
    reflect any table/view that does not have a primary key as this is a
    limitation of sqlalchemy.

    Modified from
    https://charleslavery.com/notes/sqlalchemy-reflect-tables-to-declarative.html
    """

    metadata.bind = engine

    metadata.reflect(only=should_reflect, views=True)

    define_missing_views(metadata)

    glo = globals()
    for name, obj in metadata.tables.items():
        classname = camelize_classname(name)

        try:
            glo[classname] = type(classname, (Base,), {'__table__': obj})
        except:
            logger.debug("Could not reflect table %s", name)
