import re
import random
import logging
import datetime

from sqlalchemy import desc
from sqlalchemy import Text
from sqlalchemy import Column
from sqlalchemy import String
from sqlalchemy import Integer
from sqlalchemy import DateTime
from sqlalchemy import MetaData
from sqlalchemy.dialects.mysql import LONGTEXT

from sqlalchemy.ext.automap import automap_base


logger = logging.getLogger(__name__)

metadata = MetaData()
Base = automap_base(metadata=metadata)


def camelize_classname(base, tablename, table):
    return str(tablename[0].upper() +
               re.sub(r'_(\w)', lambda m: m.group(1).upper(), tablename[1:]))


def should_reflect(tablename, *args):
    print(tablename)
    name = camelize_classname(None, tablename, None)
    return name not in globals()


def reflect(engine):
    metadata.bind = engine
    metadata.reflect(only=should_reflect)
    Base.prepare(classname_for_table=camelize_classname)
    for klass in Base.classes:
        globals()[klass.__name__] = klass


class LoopReleases(Base):
    __tablename__ = 'loop_releases'

    id = Column(String(4), primary_key=True)
    date = Column(DateTime)
    description = Column(Text)

    def __init__(self, mode='', description=''):
        self.description = description
        self.mode = mode
        self.compute_new_release_id()
        self.get_date()

    def compute_new_release_id(self, session):
        prev = session.query(LoopReleases).\
            order_by(desc(LoopReleases.date)).\
            first()
        if prev is None:
            self.id = 0.1
        elif self.mode == 'major':
            parts = prev.id.split('.')
            self.id = '.'.join([str(int(parts[0])+1), '0'])
        else:
            parts = prev.id.split('.')
            self.id = '.'.join([parts[0], str(int(parts[1])+1)])

    def get_date(self):
        self.date = datetime.datetime.now()


class MlReleases(Base):
    __tablename__ = 'ml_releases'

    id = Column(String(4), primary_key=True)
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

    id = Column(String(11), primary_key=True)
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


class NrClasses(Base):
    __tablename__ = 'nr_classes'

    id = Column(String(17), primary_key=True)
    release_id = Column(String(6),  primary_key=True)
    resolution = Column(String(4))
    handle = Column(String(5))
    version = Column(String(3))
    comment = Column(Text)

    def __init__(self, id='', increment=False, release_id='', resolution='',
                 comment=''):
        self.release_id = release_id
        self.resolution = resolution
        self.comment = comment
        if id == '':
            self.get_new_motif_id()
        elif increment is True:
            self.increment_motif_id(id)
        else:
            self.populate_fields(id)

    def generate_id(self):
        parts = [self.type, self.resolution, self.handle]
        self.id = '.'.join(['_'.join(parts), str(self.version)])

    def get_new_motif_id(self, session):
        self.type = 'NR'
        while True:
            self.handle = '%05d' % random.randrange(99999)
            handle = session.query(NrHandles).\
                filter(NrHandles.id == self.handle).\
                first()
            if handle is None:
                nr_class = session.query(NrClasses).\
                    filter(NrClasses.handle == self.handle).\
                    first()
                if nr_class is None:
                    h = NrHandles(id=self.handle)
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


class NrReleases(Base):
    __tablename__ = 'nr_releases'

    id = Column(String(6), primary_key=True)
    date = Column(DateTime)
    description = Column(Text)

    def __init__(self, mode='minor', description=''):
        self.description = description
        self.mode = mode
        self.compute_new_release_id()
        self.get_date()

    def compute_new_release_id(self, session):
        prev = session.query(NrReleases).order_by(desc(NrReleases.date)).\
            first()
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
