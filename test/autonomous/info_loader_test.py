from test import StageTest
from nose import SkipTest

from pymotifs.autonomous.info_loader import Info


class InfoLoadingTest(StageTest):
    loader_class = Info

    def setUp(self):
        super(InfoLoadingTest, self).setUp()
        self.info = self.loader.load("2AW7", "A")

    def test_loads_the_entity_id(self):
        self.assertEquals(1, self.info['entity'])

    def test_loads_database_id(self):
        self.assertEquals(4977, self.info['db_id'])

    def test_load_base_pairs(self):
        self.assertEquals(688, self.info['bp'])

    def test_loads_long_range(self):
        self.assertEquals(77, self.info['lr'])

    def test_loads_internal_cww(self):
        self.assertEquals(472, self.info['internal'])

    def test_loads_external_cww(self):
        self.assertEquals(0, self.info['external'])

    def test_loads_resolved_length_if_many_chains(self):
        self.assertEquals(35, self.loader.load('1ET4', 'A')['length'])

    def test_loads_experimental_length_if_many_chains(self):
        self.assertEquals(35, self.loader.load('1ET4', 'A')['exp_length'])

    def test_loads_resolved_length(self):
        self.assertEquals(1530, self.info['length'])

    def test_loads_experimental_length(self):
        self.assertEquals(1542, self.info['exp_length'])

    def test_loads_source_info(self):
        self.assertEquals(562, self.info['source'])

    def test_can_load_source_when_mapping_to_species(self):
        self.assertEquals(562, self.loader.load('3J01', 'A')['source'])

    def test_can_load_source_when_there_are_several(self):
        info = self.loader.load('3T4B', 'A')
        self.assertEquals(32630, info['source'])

    def test_can_load_source_when_source_is_none(self):
        self.assertEquals(None, self.loader.load('1ET4', 'A')['source'])


class MergeChainsTest(StageTest):
    loader_class = Info

    def setUp(self):
        super(MergeChainsTest, self).setUp()
        self.chains = self.loader.merge([
            self.loader.load('1ET4', 'B'),
            self.loader.load('1ET4', 'A')
        ])

    def test_can_merge_all_ids(self):
        self.assertEquals('1ET4|A,1ET4|B', self.chains['id'])

    def test_can_keep_the_pdb(self):
        self.assertEquals('1ET4', self.chains['pdb'])

    def test_can_sum_all_internal_bps(self):
        self.assertEquals(18, self.chains['summary']['internal'])

    def test_can_sum_all_bps(self):
        self.assertEquals(28, self.chains['summary']['bp'])

    def test_can_compute_external(self):
        self.assertEquals(0, self.chains['summary']['external'])

    def test_it_sums_the_resolved_length(self):
        self.assertEquals(70, self.chains['summary']['length'])

    def test_it_sums_the_experimental_length(self):
        self.assertEquals(70, self.chains['summary']['exp_length'])

    def test_it_sorts_by_length_and_name(self):
        val = [c['name'] for c in self.chains['chains']]
        self.assertEquals(['A', 'B'], val)


class MergingAutonomousChainsTest(StageTest):
    loader_class = Info

    def setUp(self):
        super(MergingAutonomousChainsTest, self).setUp()
        self.chains = self.loader.merge([
            self.loader.load('3J10', 'B'),
            self.loader.load('3J10', 'C'),
            self.loader.load('3J10', 'D')
        ])

    def test_it_sorts_the_name(self):
        self.assertEquals('3J10|D,3J10|B,3J10|C', self.chains['id'])

    def test_it_sorts_by_autonomy(self):
        val = [chain['name'] for chain in self.chains['chains']]
        self.assertEquals(['D', 'B', 'C'], val)


class CrossChainInteractionsTest(StageTest):
    loader_class = Info

    def test_gets_cross_chain_between_all_chains(self):
        raise SkipTest()
