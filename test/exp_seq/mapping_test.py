from test import StageTest
from nose import SkipTest

from pymotifs.exp_seq.mapping import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_data(self):
        raise SkipTest()

    def test_knows_if_has_no_data(self):
        raise SkipTest()

    def test_can_remove_data(self):
        raise SkipTest()


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

    def test_gets_a_correct_mapping(self):
        raise SkipTest()
