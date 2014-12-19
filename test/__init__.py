import json
import unittest as ut

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

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
