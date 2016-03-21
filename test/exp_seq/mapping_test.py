import pytest

from test import StageTest

from pymotifs.exp_seq.mapping import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_knows_if_has_data(self):
        pass

    @pytest.mark.skip()
    def test_knows_if_has_no_data(self):
        pass

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

    def setUp(self):
        super(GettingExpMappingTest, self).setUp()
        self.data = self.loader.exp_mapping('1S72', '9')

    def test_can_get_full_mapping(self):
        self.assertEquals(122, len(self.data))

    @pytest.mark.skip()
    def test_gets_a_correct_mapping(self):
        pass
