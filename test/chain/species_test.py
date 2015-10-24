from test import StageTest
from nose import SkipTest

from pymotifs.chains.species import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_has_no_data(self):
        self.assertFalse(self.loader.has_data('0GID'))

    def test_knows_has_data(self):
        self.assertTrue(self.loader.has_data('124D'))

    def test_can_remove_data(self):
        raise SkipTest()


class GettingDataTest(StageTest):
    loader_class = Loader

    def test_can_get_data_for_pdb(self):
        data = self.loader.data('1GID')
        val = [d.species_id for d in data]
        ans = [32630, 32630]
        self.assertEquals(ans, val)

    def test_assigns_synethic_when_unknown(self):
        data = self.loader.data('1ET4')
        val = [d.species_id for d in data]
        ans = [32630] * 5
        self.assertEquals(ans, val)

    def test_assigns_synthetic_when_multi(self):
        data = self.loader.data('3T4B')
        val = [entry.species_id for entry in data]
        ans = [32630]
        self.assertEquals(ans, val)

    def test_assigns_to_species(self):
        data = self.loader.data('4V6M')
        val = [e.species_id for e in data]
        ans = [562] * 5
        self.assertEquals(ans, val)

    def test_can_find_species_when_name_and_id_differ(self):
        raise SkipTest()
        # Test 4V4Q|*|DV 4L6M|*|V,4L6M|*|W
