import unittest as ut

from pymotifs.dispatcher import Dispatcher
from pymotifs.units.info import Loader
from pymotifs.units import Loader as UnitLoader


class DispatcherTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('units.info')

    def test_it_can_load_a_stage_by_name(self):
        val = self.dispatcher.get_stage('units.info')
        self.assertEquals(val, Loader)

    def test_it_can_load_an_aggregate_by_name(self):
        val = self.dispatcher.get_stage('units')
        self.assertEquals(val, UnitLoader)

    def test_it_fails_if_given_unknown_module(self):
        self.assertRaises(ImportError, self.dispatcher.get_stage, 'bob')


class LoadingTest(ut.TestCase):
    def setUp(self):
        self.dispatcher = Dispatcher('')

    def test_it_knows_if_something_is_a_loader(self):
        loader = self.dispatcher.is_loader('pymotifs.units.info')
        val = loader(Loader)
        self.assertTrue(val)

    def test_it_can_load_a_single_loader_from_many_imported(self):
        loader = self.dispatcher.is_loader('pymotifs.units')
        val = loader(UnitLoader)
        self.assertTrue(val)

    def test_it_will_skip_those_with_mismatched_name(self):
        loader = self.dispatcher.is_loader('pymotifs.units.info')
        val = loader(UnitLoader)
        self.assertFalse(val)