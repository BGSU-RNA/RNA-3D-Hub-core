from test import Session

import unittest as ut

from pymotifs.core.db import Session as SessionWrapper
from pymotifs.models import PdbInfo


class SessionWrapperTest(ut.TestCase):
    def setUp(self):
        self.session = SessionWrapper(Session)

    def tearDown(self):
        with self.session() as session:
            session.query(PdbInfo).\
                filter(PdbInfo.pdb_id.like('0%')).\
                delete(synchronize_session=False)

    def test_it_reraises_exceptions(self):
        def func():
            with self.session():
                raise AttributeError("bob")
        self.assertRaises(AttributeError, func)

    def test_it_rolls_back_on_exceptions(self):
        def func():
            with self.session() as session:
                session.add(PdbInfo(pdb_id='0000'))
                raise ValueError("Stop")

        self.assertRaises(ValueError, func)

        with self.session() as session:
            count = session.query(PdbInfo).filter_by(pdb_id='0000').count()
            self.assertEquals(0, count)

    def test_it_always_commits(self):
        with self.session() as session:
            session.add(PdbInfo(pdb_id='0000'))

        with self.session() as session:
            count = session.query(PdbInfo).filter_by(pdb_id='0000').count()
            self.assertEquals(1, count)
