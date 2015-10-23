import re
import random
import logging
import datetime

from sqlalchemy import desc
from sqlalchemy import Text
from sqlalchemy import Table
from sqlalchemy import Column
from sqlalchemy import String
from sqlalchemy import Integer
from sqlalchemy import DateTime
from sqlalchemy import MetaData
from sqlalchemy.dialects.mysql import LONGTEXT

from sqlalchemy.ext.declarative import declarative_base


logger = logging.getLogger(__name__)

metadata = MetaData()
Base = declarative_base(metadata=metadata)

PRIMARY_KEYS = {
    'correspondence_units': set(['unit_id_1', 'unit_id_2',
                                 'correspondence_id']),
    'exp_seq_pdb': set(['exp_seq_id', 'pdb_id', 'chain_id'])
}


def camelize_classname(tablename):
    return str(tablename[0].upper() +
               re.sub(r'_(\w)', lambda m: m.group(1).upper(), tablename[1:]))


def should_reflect(tablename, *args):
    name = camelize_classname(tablename)
    return name not in globals()



def reflect(engine):
    """Reflect all tables/views from the database into python. This cannot
    reflect any table/view that does not have a primary key as this is a
    limitation of sqlalchemy.

    Modified from
    https://charleslavery.com/notes/sqlalchemy-reflect-tables-to-declarative.html
    """

    metadata.bind = engine

    metadata.reflect(only=should_reflect, views=True)

    # We have to define what columns are PK for views
    Table('exp_seq_pdb', metadata,
          Column('exp_seq_id', String, primary_key=True),
          Column('pdb_id', String, primary_key=True),
          Column('chain_id', Integer, primary_key=True),
          extend_existing=True)

    glo = globals()
    for name, obj in metadata.tables.items():
        classname = camelize_classname(name)

        try:
            glo[classname] = type(classname, (Base,), {'__table__': obj})
        except:
            logger.warning("Could not reflect table %s", name)


class MlReleases(Base):
    __tablename__ = 'ml_releases'

    ml_releases_id = Column(String(4), primary_key=True)
    type = Column(String(2), primary_key=True)
    date = Column(DateTime)
    description = Column(Text)
    annotation = Column(Text)
    graphml = Column(LONGTEXT)

    def __init__(self, mode='', description='', type=''):
        self.description = description
        self.mode = mode
        self.type = type
        self.graphml = ''
        self.compute_new_release_id()
        self.get_date()

    def compute_new_release_id(self, session):
        prev = session.query(MlReleases).\
            filter(MlReleases.type == self.type).\
            order_by(desc(MlReleases.date)).\
            first()
        if prev is None:
            self.id = '0.1'
        elif self.mode == 'major':
            parts = prev.id.split('.')
            self.id = '.'.join([str(int(parts[0])+1), '0'])
        else:
            parts = prev.id.split('.')
            self.id = '.'.join([parts[0], str(int(parts[1])+1)])

    def get_date(self):
        self.date = datetime.datetime.now()


class MlMotifs(Base):
    __tablename__ = 'ml_motifs'

    ml_motifs_id = Column(String(11), primary_key=True)
    release_id = Column(String(4), primary_key=True)
    type = Column(String(2))
    handle = Column(String(5))
    version = Column(Integer)
    comment = Column(Text)

    def __init__(self, id='', release_id='', increment=False, comment='',
                 type=''):
        self.release_id = release_id
        self.comment = comment
        self.type = type

        if id == '':
            self.get_new_motif_id()
        elif increment is True:
            self.increment_motif_id(id)
        else:
            self.populate_fields(id)

    def get_new_motif_id(self, session):
        while True:
            self.handle = '%05d' % random.randrange(99999)
            motif = session.query(MlMotifs).\
                filter(MlMotifs.handle == self.handle).\
                first()
            if motif is None:
                handle = session.query(MlHandles).\
                    filter(MlHandles.id == self.handle).\
                    first()
                if handle is None:
                    h = MlHandles(id=self.handle)
                    session.add(h)
                    break
        self.version = 1
        self.id = self.type + '_' + self.handle + '.' + str(self.version)

    def increment_motif_id(self, id):
        self.handle = id[3:8]
        self.version = int(id[9:]) + 1
        self.id = self.type + '_' + self.handle + '.' + str(self.version)

    def populate_fields(self, id):
        self.id = id
        self.handle = id[3:8]
        self.version = int(id[9:])
