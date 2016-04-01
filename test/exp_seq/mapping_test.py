import pytest

from test import StageTest

from pymotifs import core
from pymotifs.exp_seq.mapping import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_data(self):
        assert self.loader.has_data('1GID') is True

    def test_knows_if_has_no_data(self):
        with pytest.raises(core.Skip):
            self.loader.has_data('0FJG')

    @pytest.mark.skip()
    def test_can_remove_data(self):
        pass


class DeterminingTheMappedChainsTest(StageTest):
    loader_class = Loader

    def test_gets_correct_chains(self):
        val = sorted(self.loader.mapped_chains('1S72'))
        self.assertEquals(['0', '9'], val)


class GettingExpMappingTest(StageTest):
    loader_class = Loader

    def test_can_get_full_mapping(self):
        val = self.loader.exp_mapping('1S72', ['9', '0'])
        self.assertEquals(2922, len(val))

    def test_can_get_full_mapping_with_duplicate_exp_seq(self):
        val = self.loader.exp_mapping('1GID', ['A', 'B'])
        self.assertEquals(158, len(val))

    def test_complains_if_cannot_get_mappings_for_all_chains(self):
        with pytest.raises(core.InvalidState):
            val = self.loader.exp_mapping('1GID', ['A', 'X'])

    @pytest.mark.skip()
    def test_gets_a_correct_mapping(self):
        pass


class DataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(DataTest, self).setUp()
        self.data = list(self.loader.data('1GID'))

    def test_it_maps_all_unit_ids(self):
        assert len(self.data) == 158 * 2
