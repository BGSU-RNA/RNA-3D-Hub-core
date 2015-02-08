from test import StageTest

from pymotifs.core import Loader
from pymotifs.models import PdbInfo


class ExampleLoader(Loader):
    def has_data(self, arg, **kwargs):
        return arg == 'known'

    def remove(self, arg):
        return True

    def data(self, arg):
        return 'bob'


class RecomputingTest(StageTest):
    loader_class = ExampleLoader

    def test_recalculates_if_given_recompute(self):
        self.assertTrue(self.loader.should_process(None, recalculate=True))

    def test_will_not_recalc_if_has_data(self):
        self.assertFalse(self.loader.should_process('known'))

    def test_will_recalc_if_missing_data(self):
        self.assertFalse(self.loader.should_process('missing'))

    def test_recalculates_if_given_on_known_data(self):
        self.assertTrue(self.loader.should_process('known', recalculate=True))

    def test_will_recalc_if_missing_data_and_recalculate(self):
        self.assertTrue(self.loader.should_process('missing',
                                                   recalculate=True))

    def test_knows_if_forced_to_recalculate(self):
        self.assertTrue(self.loader.must_recompute('missing',
                                                   recalculate=True))

    def test_knows_if_not_forced_to_recalculate(self):
        self.assertFalse(self.loader.must_recompute('missing'))

    def test_it_does_not_use_time_gap_if_not_set(self):
        self.assertFalse(self.loader.been_long_enough('missing'))


class StoringTest(StageTest):
    loader_class = ExampleLoader

    def tearDown(self):
        with self.loader.session() as session:
            session.query(PdbInfo).\
                filter(PdbInfo.id.like('0%')).\
                delete(synchronize_session=False)

    def store_and_count(self, data, **kwargs):
        self.loader.store(data, **kwargs)
        with self.loader.session() as session:
            return session.query(PdbInfo).filter(PdbInfo.id.like('0%')).count()

    def test_given_empty_list_store_nothing(self):
        self.assertEquals(0, self.store_and_count([]))

    def test_it_can_add_one_entry(self):
        self.assertEquals(1, self.store_and_count(PdbInfo(id='0000')))

    def test_it_adds_all_entries(self):
        data = [PdbInfo(id='0000'), PdbInfo(id='000A'), PdbInfo(id='000B')]
        self.assertEquals(3, self.store_and_count(data))

    def test_it_can_store_a_generator(self):
        def data():
            for id in ['0000', '000A', '000B', '000C']:
                yield PdbInfo(id=id)
        self.assertEquals(4, self.store_and_count(data()))

    def test_it_will_roll_back_if_exception_is_raised(self):
        self.assertEquals(1, self.store_and_count(PdbInfo(id='0000')))
        self.assertRaises(Exception, self.loader.store, PdbInfo(id='0000'))

    def test_it_will_not_write_if_given_dry_run(self):
        data = [PdbInfo(id='0000'), PdbInfo(id='000A'), PdbInfo(id='000B')]
        self.assertEquals(0, self.store_and_count(data, dry_run=True))

    def test_it_can_merge_if_requested(self):
        self.loader.store(PdbInfo(id='0000', resolution=10))
        with self.loader.session() as session:
                result = session.query(PdbInfo.resolution).\
                    filter_by(id='0000').first()
                self.assertEquals(10, result.resolution)

        self.loader.merge_data = True
        self.loader.store(PdbInfo(id='0000', resolution=1))

        with self.loader.session() as session:
                result = session.query(PdbInfo.resolution).\
                    filter_by(id='0000').first()
                self.assertEquals(1, result.resolution)
