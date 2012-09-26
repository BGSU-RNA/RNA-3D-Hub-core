"""

Contains declarations of NR database tables. Responsible for creation of
nr_* tables.

"""

import random
import datetime
import math
import sys

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
import sqlalchemy.exc


from get_session import get_session


"""when run from unittests, it should only use the test database"""
if 'unittest' in sys.modules:
    session, engine = get_session('test')
else:
    session, engine = get_session('dev')

session, engine = get_session('test')

Base = declarative_base()


# Utility functions
def drop_all():
    Base.metadata.drop_all(engine)

def create_all():
    Base.metadata.create_all(engine)


class NR_pdb(Base): # = loop
    """
    """
    __tablename__ = 'nr_pdbs'
    id         = Column(String(4),  primary_key=True)
    class_id   = Column(String(17), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    rep        = Column(Boolean)

    def __repr__(self):
        return "<PDB('%s','%s','%s')>" % (self.id, self.class_id, self.release_id)


class NR_class(Base): # = motif
    """
    """
    __tablename__ = 'nr_classes'

    id         = Column(String(17), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    resolution = Column(String(4))
    handle     = Column(String(5))
    version    = Column(String(3))
    comment    = Column(Text)

    def __init__(self, id='', increment=False, release_id='', resolution='', comment=''):
        self.release_id = release_id
        self.resolution = resolution
        self.comment    = comment
        if id == '':
            self.get_new_motif_id()
        elif increment is True:
            self.increment_motif_id(id)
        else:
            self.populate_fields(id)

    def generate_id(self):
        self.id = '.'.join(['_'.join([self.type,self.resolution,self.handle]),str(self.version)])

    def get_new_motif_id(self):
        self.type = 'NR'
        while True:
            self.handle = '%05d' % random.randrange(99999)
            if session.query(NR_handle).filter(NR_handle.id==self.handle).first() is None:
                if session.query(NR_class).filter(NR_class.handle==self.handle).first() is None:
                    h = NR_handle(id=self.handle)
                    session.add(h)
                    break
        self.version = 1
        self.generate_id()

    def populate_fields(self, id):
        self.id = id
        [self.type, self.resolution, handle_version] = id.split('_')
        [self.handle, version] = handle_version.split('.')
        self.version = int(version)

    def increment_motif_id(self, id):
        self.populate_fields(id)
        self.version = int(self.version) + 1
        self.generate_id()

    def __repr__(self):
        return "<NR_class('%s','%s','%s')>" % (self.id, self.type, self.release_id)


class NR_handle(Base):
    """
    """
    __tablename__ = 'nr_handles'
    id = Column(String(5), primary_key=True)


class NR_release(Base):
    """
    mode = minor/major/reuse. Release mode
    """
    __tablename__ = 'nr_releases'

    id          = Column(String(4), primary_key=True)
    date        = Column(DateTime)
    description = Column(Text)

    def __init__(self, mode='minor', description=''):
        self.description = description
        self.mode = mode
        self.compute_new_release_id()
        self.get_date()

    def __repr__(self):
        return "<NR_release('%s','%s','%s')>" % (self.id, self.date, self.description)

    def compute_new_release_id(self):
        prev = session.query(NR_release).order_by(desc(NR_release.date)).first()
        if prev is None:
            self.id = '0.1'
        elif self.mode == 'major':
            parts = prev.id.split('.')
            self.id = '.'.join([str(int(parts[0])+1), '0'])
        elif self.mode == 'reuse':
            self.id = prev.id
        else:
            parts = prev.id.split('.')
            self.id = '.'.join([parts[0], str(int(parts[1])+1)])

    def get_date(self):
        self.date = datetime.datetime.now()


class NR_setdiff(Base):
    """
    """
    __tablename__ = 'nr_set_diff'

    nr_class1     = Column(String(17), primary_key=True)
    nr_class2     = Column(String(17), primary_key=True)
    release_id    = Column(String(4),  primary_key=True)
    intersection  = Column(Text)
    overlap       = Column(Float)
    one_minus_two = Column(Text)
    two_minus_one = Column(Text)

    def __repr__(self):
        return "<NRSetDiff('%s','%s','%s')>" % (self.nr_class1, self.nr_class2, self.release_id)


class NR_parents(Base):
    """
    """
    __tablename__ = 'nr_parents'

    class_id   = Column(String(17), primary_key=True)
    release_id = Column(String(4),  primary_key=True)
    parents    = Column(Text)

    def __repr__(self):
        return "<NR_parents('%s','%s','%s')>" % (self.motif_id, self.release_id, self.parents)


class NR_release_diff(Base):
    """
    """
    __tablename__ = 'nr_release_diff'

    nr_release_id1 = Column(String(4), primary_key=True)
    nr_release_id2 = Column(String(4), primary_key=True)
    resolution     = Column(String(4), primary_key=True)
    direct_parent  = Column(Boolean)
    added_groups   = Column(Text)
    removed_groups = Column(Text)
    updated_groups = Column(Text)
    same_groups    = Column(Text)
    added_pdbs     = Column(Text)
    removed_pdbs   = Column(Text)
    num_added_groups   = Column(Integer)
    num_removed_groups = Column(Integer)
    num_updated_groups = Column(Integer)
    num_same_groups    = Column(Integer)
    num_added_pdbs     = Column(Integer)
    num_removed_pdbs   = Column(Integer)

    def __repr__(self):
        return "<NRReleaseDiff('%s','%s')>" % (self.nr_release_id1, self.nr_release_id2)


Base.metadata.create_all(engine)
