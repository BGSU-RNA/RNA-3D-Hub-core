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
        val = sorted([(entry.chain_id, entry.species_id) for entry in data],
                     key=lambda e: e[0])
        ans = [(492, 32630), (493, 32630)]
        self.assertEquals(ans, val)

    def test_assigns_synethic_when_unknown(self):
        data = self.loader.data('1ET4')
        val = sorted([(entry.chain_id, entry.species_id) for entry in data],
                     key=lambda e: e[0])
        ans = [(301, 32630), (302, 32630), (303, 32630), (304, 32630),
               (305, 32630)]
        self.assertEquals(ans, val)

    def test_assigns_synthetic_when_multi(self):
        data = self.loader.data('3T4B')
        val = [(entry.chain_id, entry.species_id) for entry in data]
        ans = [(18248, 32630)]
        self.assertEquals(ans, val)

    def test_assigns_to_species(self):
        data = self.loader.data('3J01')
        val = sorted([(entry.chain_id, entry.species_id) for entry in data],
                     key=lambda e: e[0])
        ans = [(14398, 562), (14399, 562)]
        self.assertEquals(ans, val)
