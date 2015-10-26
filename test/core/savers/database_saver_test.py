from test import CONFIG
from test import StageTest
from test import Session as SessionMaker

from pymotifs.core.db import Session
from pymotifs.core.exceptions import InvalidState
from pymotifs.core.savers import DatabaseSaver
from pymotifs.models import PdbInfo

Session = Session(SessionMaker)


class SimpleDatabaseSavingTest(StageTest):
    def setUp(self):
        self.saver = DatabaseSaver(CONFIG, Session)

    def tearDown(self):
        with Session() as session:
            session.query(PdbInfo).\
                filter(PdbInfo.pdb_id.like('0%')).\
                delete(synchronize_session=False)

    def store_and_count(self, name, data, **kwargs):
        self.saver(name, data, **kwargs)
        return self.count()

    def count(self):
        with Session() as session:
            return session.query(PdbInfo).\
                filter(PdbInfo.pdb_id.like('0%')).\
                count()

    def test_it_will_complain_given_nothing(self):
        self.assertRaises(InvalidState, self.saver, '0000', [])

    def test_it_can_add_one_entry(self):
        val = self.store_and_count('0000', PdbInfo(pdb_id='0000'))
        self.assertEquals(1, val)

    def test_it_can_save_one_dict(self):
        self.saver.table = PdbInfo
        data = {'pdb_id': '0000', 'resolution': 10}
        self.assertEquals(1, self.store_and_count('0000', data))

    def test_it_can_save_several_dicts(self):
        self.saver.table = PdbInfo
        data = [{'pdb_id': '0000', 'resolution': 10},
                {'pdb_id': '000A', 'resolution': 2}]
        self.assertEquals(2, self.store_and_count('0000', data))

    def test_it_adds_all_entries(self):
        data = [PdbInfo(pdb_id='0000'), PdbInfo(pdb_id='000A'),
                PdbInfo(pdb_id='000B')]
        self.assertEquals(3, self.store_and_count('0000', data))

    def test_it_can_store_a_generator(self):
        def data():
            for id in ['0000', '000A', '000B', '000C']:
                yield PdbInfo(pdb_id=id)
        self.assertEquals(4, self.store_and_count('0000', data()))

    def test_it_will_complain_if_exception_in_generator(self):
        def data():
            for index, id in enumerate(['0000', '000A', '000B', '000C']):
                if index == 1:
                    raise ValueError()
                yield PdbInfo(pdb_id=id)
        saver = DatabaseSaver(CONFIG, Session)
        self.assertRaises(ValueError, saver, '0000', data())

    def test_it_will_not_write_if_given_dry_run(self):
        data = [PdbInfo(pdb_id='0000'), PdbInfo(pdb_id='000A'),
                PdbInfo(pdb_id='000B')]
        self.assertEquals(0, self.store_and_count('0000', data, dry_run=True))

    def test_it_can_merge_if_requested(self):
        saver = DatabaseSaver(CONFIG, Session)
        saver.merge = True
        saver('0000', PdbInfo(pdb_id='0000', resolution=10))
        saver('0000', PdbInfo(pdb_id='0000', resolution=9))

        with self.saver.session() as session:
                result = session.query(PdbInfo.resolution).\
                    filter_by(pdb_id='0000').one()
                self.assertEquals(9, result.resolution)
