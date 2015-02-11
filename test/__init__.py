import os
import json
import unittest as ut
from functools import wraps

from nose import SkipTest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from fr3d.cif.reader import Cif

import pymotifs.models as models


with open('conf/test.json', 'rb') as raw:
    config = json.load(raw)

engine = create_engine(config['db']['uri'])
models.reflect(engine)
Session = sessionmaker(bind=engine)


def which(program):

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def skip_without_matlab(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        if not which('matlab'):
            raise SkipTest("Skipping without matlab")
        return func(*args, **kwargs)
    return wrapper


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
        else:
            cls.cif = None
            cls.structure = None

    def setUp(self):
        super(CifStageTest, self).setUp
        self.structure = self.__class__.structure
        self.cif = self.__class__.cif
