import pytest

from test import StageTest

from pymotifs import core
from pymotifs.exp_seq.mapping import Loader
from pymotifs.exp_seq.mapping import MappedChain


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_data(self):
        assert self.loader.has_data('1GID') is True

    def test_knows_if_has_no_data(self):
        assert self.loader.has_data('0FJG') is False

    @pytest.mark.skip()
    def test_can_remove_data(self):
        pass


class DeterminingTheMappedChainsTest(StageTest):
    loader_class = Loader

    def test_gets_correct_chains(self):
        val = sorted(self.loader.mapped_chains('1S72'))
        assert len(val) == 2
        assert sorted(v.name for v in val) == ['0', '9']


class GettingExpMappingTest(StageTest):
    loader_class = Loader

    def exp_mapping(self, pdb, *chains):
        mapped = [MappedChain(id=None, chain_id=None, name=c) for c in chains]
        return self.loader.exp_mapping(pdb, mapped)

    def test_can_get_full_mapping(self):
        val = self.exp_mapping('1S72', '0', '9')
        self.assertEquals(2922 + 122, len(val))

    def test_can_get_full_mapping_with_duplicate_exp_seq(self):
        val = self.exp_mapping('1GID', 'A', 'B')
        self.assertEquals(158 * 2, len(val))

    def test_complains_if_cannot_get_mappings_for_all_chains(self):
        with pytest.raises(core.InvalidState):
            self.exp_mapping('1GID', 'A', 'X')

    def test_gets_a_correct_mapping(self):
        val = self.exp_mapping('1GID', 'A', 'B')
        ans = []
        for chain in ['A', 'B']:
            for index in range(158):
                ans.append((chain, index))
        assert ans == sorted(val.keys())


class DataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(DataTest, self).setUp()
        self.data = list(self.loader.data('1GID'))

    def test_it_maps_all_unit_ids(self):
        assert len(self.data) == 158 * 2
