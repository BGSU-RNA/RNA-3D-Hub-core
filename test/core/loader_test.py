from test import StageTest

from pymotifs.core import Loader


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
