from test import StageTest
from nose import SkipTest

from pymotifs.nr.groups import Grouper


class InfoLoadingTest(StageTest):
    loader_class = Grouper

    def setUp(self):
        super(InfoLoadingTest, self).setUp()
        self.info = self.loader.info("2AW7", "A")

    def test_load_base_pairs(self):
        raise SkipTest()
        self.assertEquals(689, self.info['bp'])

    def test_loads_long_range(self):
        raise SkipTest()
        self.assertEquals(None, self.info['lr'])

    def test_loads_internal_cww(self):
        raise SkipTest()
        self.assertEquals(None, self.info['internal'])

    def test_loads_external_cww(self):
        self.assertEquals(0, self.info['external'])

    def test_loads_database_id(self):
        self.assertEquals(4977, self.info['db_id'])

    def test_loads_resolved_length(self):
        raise SkipTest()
        self.assertEquals(None, self.info['length'])

    def test_loads_experimental_length(self):
        raise SkipTest()
        self.assertEquals(None, self.info['length'])

    def test_loads_source_info(self):
        raise SkipTest()
        self.assertEquals(None, self.info['source'])


class MergeChainsTest(StageTest):
    loader_class = Grouper

    def test_can_merge_all_ids(self):
        raise SkipTest()

    def test_can_keep_the_pdb(self):
        raise SkipTest()

    def test_can_merge_all_chains(self):
        raise SkipTest()

    def test_can_sum_all_internal(self):
        raise SkipTest()

    def test_can_sum_all_bps(self):
        raise SkipTest()

    def test_can_compute_external(self):
        raise SkipTest()


class ChainsTest(StageTest):
    loader_class = Grouper

    def test_can_find_simple_groups(self):
        val = self.loader.chains("2AW7")
        self.assertEquals(1, len(val))

    def test_can_group_several_chains(self):
        val = self.loader.chains('1ET4')
        self.assertEquals(3, len(val))

    def test_it_can_create_proper_mappings(self):
        grouping = self.loader.chains('1ET4')
        val = sorted([sorted(group['name']) for group in grouping])
        ans = [['A', 'B'], ['C', 'D'], ['E']]
        self.assertEquals(ans, val)
