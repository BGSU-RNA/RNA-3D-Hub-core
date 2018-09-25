from test import StageTest

from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod


class Simple(core.SimpleLoader):
    def query(self, session, pdb):
        return session.query(mod.UnitInfo).\
            join(mod.PdbInfo, mod.PdbInfo.pdb_id == mod.UnitInfo.pdb_id).\
            filter(mod.PdbInfo.pdb_id == pdb)

    def data(self, pdb, **kwargs):
        return [pdb]


class HasDataTest(StageTest):
    loader_class = Simple

    def test_it_can_tell_if_data_exists(self):
        self.assertTrue(self.loader.has_data('1GID'))

    def test_it_can_tell_if_no_data_exists(self):
        self.assertFalse(self.loader.has_data('0GID'))


class RemovingTest(StageTest):
    loader_class = Simple

    def setUp(self):
        super(RemovingTest, self).setUp()
        with self.loader.session() as session:
            session.add(mod.PdbInfo(pdb_id='0000'))

    def tearDown(self):
        with self.loader.session() as session:
            session.query(mod.PdbInfo).\
                filter_by(pdb_id='0000').\
                delete(synchronize_session=False)

    def count(self):
        with self.loader.session() as session:
            query = session.query(mod.UnitInfo).\
                filter(mod.UnitInfo.pdb_id == '0000')
            return query.count()

    def dummy(self, count):

        def generator():
            pattern = '0000|1|C||%s'
            for index in xrange(count):
                yield pattern % index

        for chunk in ut.grouper(1000, generator()):
            with self.loader.session() as session:
                for d in chunk:
                    session.add(mod.UnitInfo(unit_id=d, pdb_id='0000'))

        self.assertEquals(count, self.count())

    def test_it_can_remove_an_entry(self):
        self.dummy(1)
        self.loader.remove('0000')
        self.assertEquals(0, self.count())

    def test_it_can_remove_many_entries(self):
        self.dummy(10000)
        self.loader.remove('0000')
        self.assertEquals(0, self.count())

    def test_it_can_handle_nothing_to_remove(self):
        self.assertTrue(self.loader.remove('0000'))

    def test_it_does_nothing_given_dry_run(self):
        self.dummy(1)
        self.loader.remove('0000', dry_run=True)
        self.assertEquals(1, self.count())
