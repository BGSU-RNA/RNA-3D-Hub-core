from nose import SkipTest

from test import StageTest

from pymotifs.nr.groups import Grouper


class InfoLoadingTest(StageTest):
    loader_class = Grouper

    def setUp(self):
        super(InfoLoadingTest, self).setUp()
        self.info = self.loader.info("2AW7", "A")

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
        self.assertEquals(35, self.loader.info('1ET4', 'A')['length'])

    def test_loads_experimental_length_if_many_chains(self):
        self.assertEquals(35, self.loader.info('1ET4', 'A')['exp_length'])

    def test_loads_resolved_length(self):
        self.assertEquals(1530, self.info['length'])

    def test_loads_experimental_length(self):
        self.assertEquals(1542, self.info['exp_length'])

    def test_loads_source_info(self):
        raise SkipTest()
        self.assertEquals([562], self.info['source'])

    def test_can_load_source_when_mapping_to_species(self):
        raise SkipTest()
        self.assertEquals([562], self.loader.info('3J01', 'A')['source'])

    def test_can_load_source_when_there_are_several(self):
        raise SkipTest()
        info = self.loader.info('3T4B', 'A')
        self.assertEquals([32630, 11103], info['source'])

    def test_can_load_source_when_source_is_none(self):
        raise SkipTest()
        self.assertEquals(None, self.loader.info('1ET4', 'A')['source'])


class MergeChainsTest(StageTest):
    loader_class = Grouper

    def setUp(self):
        super(MergeChainsTest, self).setUp()
        self.chains = self.loader.merge_chains([
            self.loader.info('1ET4', 'B'),
            self.loader.info('1ET4', 'A')
        ])

    def test_can_merge_all_ids(self):
        self.assertEquals('1ET4|A,1ET4|B', self.chains['id'])

    def test_can_keep_the_pdb(self):
        self.assertEquals('1ET4', self.chains['pdb'])

    def test_can_merge_all_chain_names(self):
        self.assertEquals(['A', 'B'], self.chains['name'])

    def test_merges_db_id(self):
        self.assertEquals([301, 302], self.chains['db_id'])

    def test_merges_the_entity_id(self):
        self.assertEquals([1, 1], self.chains['entity'])

    def test_will_sort_entity_ids(self):
        chains = self.loader.merge_chains([
            self.loader.info('1FEU', 'A'),
            self.loader.info('1FEU', 'B')
        ])
        self.assertEquals([1, 3], chains['entity'])

    def test_can_sum_all_internal_bps(self):
        self.assertEquals(18, self.chains['internal'])

    def test_can_sum_all_bps(self):
        self.assertEquals(28, self.chains['bp'])

    def test_can_compute_external(self):
        self.assertEquals(0, self.chains['external'])

    def test_it_sums_the_resolved_length(self):
        self.assertEquals(70, self.chains['length'])

    def test_it_sums_the_experimental_length(self):
        self.assertEquals(70, self.chains['exp_length'])


class BestChainsTest(StageTest):
    loader_class = Grouper

    def test_selects_based_upon_bp_per_nt_first(self):
        chains = [
            {'bp': 12, 'length': 4, 'name': 'C'},
            {'bp': 13, 'length': 4, 'name': 'A'},
        ]
        val = self.loader.best(chains)['name']
        self.assertEquals('A', val)

    def test_selects_based_upon_length_second(self):
        chains = [
            {'bp': 0, 'length': 10, 'name': 'C'},
            {'bp': 0, 'length': 8, 'name': 'A'},
        ]
        val = self.loader.best(chains)['name']
        self.assertEquals('C', val)

    def test_selects_based_upon_name_last(self):
        chains = [
            {'bp': 10, 'length': 4, 'name': 'C'},
            {'bp': 10, 'length': 4, 'name': 'A'},
        ]
        val = self.loader.best(chains)['name']
        self.assertEquals('A', val)

    def test_works_even_if_there_are_no_bps(self):
        chains = [
            {'bp': 0, 'length': 4, 'name': 'C'},
            {'bp': 0, 'length': 4, 'name': 'A'},
        ]
        val = self.loader.best(chains)['name']
        self.assertEquals('A', val)

    def test_works_with_lists_of_names(self):
        chains = [
            {'bp': 10, 'length': 4, 'name': ['C', 'D']},
            {'bp': 10, 'length': 4, 'name': ['A', 'D']},
        ]
        val = self.loader.best(chains)['name']
        self.assertEquals(['A', 'D'], val)


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
        val = sorted([group['name'] for group in grouping])
        ans = [['A', 'B'], ['C', 'D'], ['E']]
        self.assertEquals(ans, val)

    def test_it_can_process_1F5H(self):
        val = self.loader.chains('1F5H')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1F5H(self):
        grouping = self.loader.chains('1F5H')
        val = sorted([group['name'] for group in grouping])
        ans = [['A', 'B']]
        self.assertEquals(ans, val)

    def test_it_can_process_pseudoknotted_chains(self):
        val = self.loader.chains('1F5U')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1F5U(self):
        grouping = self.loader.chains('1F5U')
        val = sorted([group['name'] for group in grouping])
        ans = [['A', 'B']]
        self.assertEquals(ans, val)

    def test_it_can_process_paired_strands(self):
        val = self.loader.chains('1EKD')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1EKD(self):
        grouping = self.loader.chains('1EKD')
        val = sorted([group['name'] for group in grouping])
        ans = [['A', 'B']]
        self.assertEquals(ans, val)

    def test_it_can_process_identical_chains(self):
        val = self.loader.chains('1EIY')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1EIY(self):
        grouping = self.loader.chains('1EIY')
        val = sorted([group['name'] for group in grouping])
        ans = [['C']]
        self.assertEquals(ans, val)

    def test_it_can_process_paired_identical_chains(self):
        val = self.loader.chains('1FEU')
        self.assertEquals(2, len(val))

    def test_it_can_select_the_correct_chains_of_1FEU(self):
        grouping = self.loader.chains('1FEU')
        val = sorted([group['name'] for group in grouping])
        ans = [['B', 'C'], ['E', 'F']]
        self.assertEquals(ans, val)

    def test_can_process_a_messy_nmr(self):
        val = self.loader.chains('1FCW')
        self.assertEquals(4, len(val))

    def test_can_select_the_correct_chains_of_1FCW(self):
        grouping = self.loader.chains('1FCW')
        val = sorted([group['name'] for group in grouping])
        ans = [['A', 'B'], ['C'],  ['D'], ['E']]
        self.assertEquals(ans, val)


class NrChainsTest(StageTest):
    loader_class = Grouper

    def test_can_select_nr_chains_with_one_chain(self):
        val = self.loader.nr_chains("2AW7")
        self.assertEquals(1, len(val))

    def test_can_select_nr_chains_with_duplicates(self):
        val = self.loader.nr_chains('1EIY')
        self.assertEquals(1, len(val))

    def test_it_can_process_paired_identical_chains(self):
        val = self.loader.nr_chains('1FEU')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1FEU(self):
        grouping = self.loader.nr_chains('1FEU')
        val = sorted([group['name'] for group in grouping])
        ans = [['B', 'C']]
        self.assertEquals(ans, val)

    def test_can_process_a_messy_nmr(self):
        val = self.loader.nr_chains('1FCW')
        self.assertEquals(2, len(val))

    def test_can_process_several_distinct_chains(self):
        val = self.loader.nr_chains('1ET4')
        self.assertEquals(2, len(val))

    def test_it_selects_correct_chains(self):
        grouping = self.loader.nr_chains('1ET4')
        val = sorted([group['name'] for group in grouping])
        ans = [['A', 'B'], ['E']]
        self.assertEquals(ans, val)

    def test_it_can_process_1F5H(self):
        val = self.loader.nr_chains('1F5H')
        self.assertEquals(1, len(val))

    def test_it_can_process_pseudoknotted_chains(self):
        val = self.loader.nr_chains('1F5U')
        self.assertEquals(1, len(val))
