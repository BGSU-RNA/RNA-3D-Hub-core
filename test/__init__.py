import json
import unittest as ut

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from fr3d.cif.reader import Cif

import pymotifs.models as models

with open('conf/test.json', 'rb') as raw:
    config = json.load(raw)

engine = create_engine(config['db']['uri'])
models.reflect(engine)
Session = sessionmaker(bind=engine)


class StageTest(ut.TestCase):
    loader_class = None

    def setUp(self):
        if self.loader_class:
            self.loader = self.loader_class(config, Session)


class QueryUtilTest(ut.TestCase):
    query_class = None

    def setUp(self):
        if self.query_class:
            self.db_obj = self.query_class(Session)


class CifStageTest(StageTest):
    filename = None

    @classmethod
    def setUpClass(cls):
        if cls.filename:
            with open(cls.filename, 'rb') as raw:
                cls.cif = Cif(raw)
                cls.structure = cls.cif.structure()

    def setUp(self):
        super(CifStageTest, self).setUp
        self.structure = self.__class__.structure
        self.cif = self.__class__.cif
