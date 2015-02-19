from test import StageTest
from nose import SkipTest

from pymotifs.chains.species import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_has_no_data(self):
        self.assertFalse(self.loader.has_data('0GID'))

    def test_knows_has_data(self):
        raise SkipTest()

    def test_can_remove_data(self):
        raise SkipTest()


class GettingDataTest(StageTest):
    loader_class = Loader

    def test_can_get_data_for_pdb(self):
        data = self.loader.data('1GID')
        val = [(entry.id, entry.species_id) for entry in data]
        ans = [('A', 32630), ('B', 32630)]
        self.assertEquals(ans, val)
