from test import StageTest
from test import skip_without_matlab
from nose import SkipTest

from pymotifs.ife.info_loader import Info


class InfoLoadingTest(StageTest):
    loader_class = Info

    def setUp(self):
        super(InfoLoadingTest, self).setUp()
        self.info = self.loader.load("4V4Q", "AA")

    def test_loads_the_entity_id(self):
        self.assertEquals(1, self.info['entity'])

    def test_loads_database_id(self):
        self.assertTrue(self.info['db_id'])

    @skip_without_matlab
    def test_load_base_pairs(self):
        self.assertEquals(680, self.info['bp'])

    @skip_without_matlab
    def test_loads_long_range(self):
        self.assertEquals(71, self.info['lr'])

    @skip_without_matlab
    def test_loads_internal_cww(self):
        self.assertEquals(472, self.info['internal'])

    @skip_without_matlab
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
        self.assertEquals(562, self.loader.load('4V6M', 'AA')['source'])

    def test_can_load_source_when_there_are_several(self):
        info = self.loader.load('3T4B', 'A')
        self.assertEquals(None, info['source'])

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
        self.assertEquals('1ET4|A+1ET4|B', self.chains['id'])

    def test_can_keep_the_pdb(self):
        self.assertEquals('1ET4', self.chains['pdb'])

    @skip_without_matlab
    def test_can_sum_all_internal_bps(self):
        self.assertEquals(18, self.chains['summary']['internal'])

    @skip_without_matlab
    def test_can_sum_all_bps(self):
        self.assertEquals(28, self.chains['summary']['bp'])

    @skip_without_matlab
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
            self.loader.load('4V6R', 'AB'),
            self.loader.load('4V6R', 'AC'),
            self.loader.load('4V6R', 'AD')
        ])

    def test_it_sorts_the_name(self):
        self.assertEquals('4V6R|AD+4V6R|AB+4V6R|AC', self.chains['id'])

    def test_it_sorts_by_autonomy(self):
        val = [chain['name'] for chain in self.chains['chains']]
        self.assertEquals(['AD', 'AB', 'AC'], val)


class CrossChainInteractionsTest(StageTest):
    loader_class = Info

    def test_gets_cross_chain_between_all_chains(self):
        raise SkipTest()


class LoadingIssuesTest(StageTest):
    loader_class = Info

    def test_loads_correct_chains_for_1G59(self):
        val = self.loader.rna_chains('1G59')
        self.assertEquals(['B', 'D'], val)
