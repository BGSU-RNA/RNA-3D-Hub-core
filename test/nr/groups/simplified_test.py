from test import StageTest
from nose import SkipTest

from pymotifs.nr.groups.simplified import Grouper
from pymotifs.nr.groups.simplified import ranking_key


class RankingChainsTest(StageTest):
    loader_class = Grouper

    def sort(self, chains):
        ordered = sorted(chains, key=ranking_key)
        return [c['name'] for c in ordered]

    def best(self, chains):
        return self.sort(chains)[0]

    def test_selects_based_upon_bp_per_nt_first(self):
        chains = [
            {'bp': 12, 'length': 4, 'name': 'C'},
            {'bp': 13, 'length': 4, 'name': 'A'},
        ]
        val = self.best(chains)
        self.assertEquals('A', val)

    def test_selects_based_upon_length_second(self):
        chains = [
            {'bp': 0, 'length': 10, 'name': 'C'},
            {'bp': 0, 'length': 8, 'name': 'A'},
        ]
        val = self.best(chains)
        self.assertEquals('C', val)

    def test_perfers_some_with_bps_over_missing(self):
        chains = [
            {'bp': 0, 'length': 10, 'name': 'C'},
            {'bp': 1, 'length': 8, 'name': 'A'},
        ]
        val = self.best(chains)
        self.assertEquals('A', val)

    def test_selects_based_upon_name_last(self):
        chains = [
            {'bp': 10, 'length': 4, 'name': 'C'},
            {'bp': 10, 'length': 4, 'name': 'A'},
        ]
        val = self.best(chains)
        self.assertEquals('A', val)

    def test_works_even_if_there_are_no_bps(self):
        chains = [
            {'bp': 0, 'length': 4, 'name': 'C'},
            {'bp': 0, 'length': 4, 'name': 'A'},
        ]
        val = self.best(chains)
        self.assertEquals('A', val)

    def test_works_with_lists_of_names(self):
        chains = [
            {'bp': 10, 'length': 4, 'name': ['C', 'D']},
            {'bp': 10, 'length': 4, 'name': ['A', 'D']},
        ]
        val = self.best(chains)
        self.assertEquals(['A', 'D'], val)


class AutonomousChainsTest(StageTest):
    loader_class = Grouper

    def grouped_names(self, pdb):
        names = []
        for group in self.loader.chains(pdb):
            names.append([chain['name'] for chain in group['chains']])
        return names

    def test_can_find_simple_groups(self):
        val = self.loader.chains("124D")
        self.assertEquals(1, len(val))

    def test_can_group_several_chains(self):
        val = self.loader.chains('1ET4')
        self.assertEquals(5, len(val))

    def test_it_can_create_proper_mappings(self):
        val = self.grouped_names('1ET4')
        ans = [['A'], ['B'], ['C'], ['D'], ['E']]
        self.assertEquals(ans, val)

    def test_it_can_process_1F5H(self):
        val = self.loader.chains('1F5H')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1F5H(self):
        val = self.grouped_names('1F5H')
        ans = [['A', 'B']]
        self.assertEquals(ans, val)

    def test_it_can_process_pseudoknotted_chains(self):
        val = self.loader.chains('1F5U')
        self.assertEquals(2, len(val))

    def test_it_can_select_the_correct_chains_of_1F5U(self):
        val = self.grouped_names('1F5U')
        ans = [['A'], ['B']]
        self.assertEquals(ans, val)

    def test_it_can_process_paired_strands(self):
        val = self.loader.chains('1EKD')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1EKD(self):
        val = self.grouped_names('1EKD')
        ans = [['A', 'B']]
        self.assertEquals(ans, val)

    def test_it_can_process_identical_chains(self):
        val = self.loader.chains('1EIY')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1EIY(self):
        val = self.grouped_names('1EIY')
        ans = [['C']]
        self.assertEquals(ans, val)

    def test_it_can_process_paired_identical_chains(self):
        val = self.loader.chains('1FEU')
        self.assertEquals(2, len(val))

    def test_it_can_select_the_correct_chains_of_1FEU(self):
        val = self.grouped_names('1FEU')
        ans = [['B', 'C'], ['E', 'F']]
        self.assertEquals(ans, val)

    def test_can_process_a_messy_nmr(self):
        val = self.loader.chains('1FCW')
        self.assertEquals(5, len(val))

    def test_can_select_the_correct_chains_of_1FCW(self):
        val = self.grouped_names('1FCW')
        ans = [['A'], ['B'], ['C'],  ['D'], ['E']]
        self.assertEquals(ans, val)

    def test_can_select_nr_chains_with_one_chain(self):
        val = self.loader.chains("124D")
        self.assertEquals(1, len(val))

    def test_can_select_nr_chains_with_duplicates(self):
        val = self.loader.chains('1EIY')
        self.assertEquals(1, len(val))

    def test_can_process_several_distinct_chains(self):
        val = self.loader.chains('1ET4')
        self.assertEquals(5, len(val))

#     def test_it_selects_correct_chains(self):
#         grouping = self.loader.chains('1ET4')
#         val = sorted([group['name'] for group in grouping])
#         ans = [['A', 'B'], ['E']]
        # self.assertEquals(ans, val)

    def test_can_process_simple_structure(self):
        val = self.grouped_names('1G59')
        ans = [['B'], ['D']]
        self.assertEquals(ans, val)

    def test_it_can_process_a_triplet(self):
        val = self.grouped_names('4V42')
        ans = [['A1', 'AB', 'AC'], ['BA'], ['BD']]
        self.assertEquals(ans, val)


class SimilarSpeciesTest(StageTest):
    loader_class = Grouper

    def test_says_different_if_disagree(self):
        raise SkipTest()

    def test_ignores_difference_if_synethic(self):
        raise SkipTest()

    def test_ignores_difference_if_not_set(self):
        raise SkipTest()

    def test_uses_only_difference_in_first_chain(self):
        raise SkipTest()
