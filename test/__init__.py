import os
import unittest as ut

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from fr3d.cif.reader import Cif

import pymotifs.models as models
from pymotifs.config import load as config_loader

from pymotifs.utils.matlab import exists as has_matlab

CONFIG = config_loader('conf/test.json', )
engine = create_engine(CONFIG['db']['uri'])
models.reflect(engine)
Session = sessionmaker(bind=engine)


skip_without_matlab = pytest.mark.skipif(
    has_matlab() is True,
    reason="No matlab installed")


class StageTest(ut.TestCase):
    loader_class = None

    def setUp(self):
        os.chdir(CONFIG['locations']['base'])
        if self.loader_class:
            self.loader = self.loader_class(CONFIG, Session)


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
        else:
            cls.cif = None
            cls.structure = None

    def setUp(self):
        super(CifStageTest, self).setUp()
        self.structure = self.__class__.structure
        self.cif = self.__class__.cif
