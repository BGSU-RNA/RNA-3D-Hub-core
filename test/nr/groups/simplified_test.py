import functools as ft

import pytest

from test import StageTest

from pymotifs.nr.groups.simplified import Grouper
from pymotifs.nr.groups.simplified import ranking_key


def basic_data(id, length=100, resolution=2.0, method='X-RAY DIFFRACTION',
               species=512):
    return {
        'id': str(id),
        'db_id': id,
        'chains': [{'db_id': id}],
        'length': length,
        'resolution': resolution,
        'method': method,
        'species': species,
    }


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
        for group in self.loader.ifes(pdb):
            names.append([chain['name'] for chain in group['chains']])
        return names

    def test_can_find_simple_groups(self):
        val = self.loader.ifes("124D")
        self.assertEquals(1, len(val))

    def test_can_group_several_chains(self):
        val = self.loader.ifes('1ET4')
        self.assertEquals(5, len(val))

    def test_it_can_create_proper_mappings(self):
        val = self.grouped_names('1ET4')
        ans = [['A'], ['B'], ['C'], ['D'], ['E']]
        self.assertEquals(ans, val)

    def test_it_can_process_1F5H(self):
        val = self.loader.ifes('1F5H')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1F5H(self):
        val = self.grouped_names('1F5H')
        ans = [['A', 'B']]
        self.assertEquals(ans, val)

    def test_it_can_process_pseudoknotted_chains(self):
        val = self.loader.ifes('1F5U')
        self.assertEquals(2, len(val))

    def test_it_can_select_the_correct_chains_of_1F5U(self):
        val = self.grouped_names('1F5U')
        ans = [['A'], ['B']]
        self.assertEquals(ans, val)

    def test_it_can_process_paired_strands(self):
        val = self.loader.ifes('1EKD')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1EKD(self):
        val = self.grouped_names('1EKD')
        ans = [['A', 'B']]
        self.assertEquals(ans, val)

    def test_it_can_process_identical_chains(self):
        val = self.loader.ifes('1EIY')
        self.assertEquals(1, len(val))

    def test_it_can_select_the_correct_chains_of_1EIY(self):
        val = self.grouped_names('1EIY')
        ans = [['C']]
        self.assertEquals(ans, val)

    def test_it_can_process_paired_identical_chains(self):
        val = self.loader.ifes('1FEU')
        self.assertEquals(2, len(val))

    def test_it_can_select_the_correct_chains_of_1FEU(self):
        val = self.grouped_names('1FEU')
        ans = [['C', 'B'], ['F', 'E']]
        self.assertEquals(ans, val)

    def test_can_process_a_messy_nmr(self):
        val = self.loader.ifes('1FCW')
        self.assertEquals(5, len(val))

    def test_can_select_the_correct_chains_of_1FCW(self):
        val = self.grouped_names('1FCW')
        ans = [['A'], ['B'], ['C'],  ['D'], ['E']]
        self.assertEquals(ans, val)

    def test_can_select_nr_chains_with_one_chain(self):
        val = self.loader.ifes("124D")
        self.assertEquals(1, len(val))

    def test_can_select_nr_chains_with_duplicates(self):
        val = self.loader.ifes('1EIY')
        self.assertEquals(1, len(val))

    def test_can_process_several_distinct_chains(self):
        val = self.loader.ifes('1ET4')
        self.assertEquals(5, len(val))

    def test_does_not_join_by_non_cww(self):
        grouping = self.loader.ifes('1ET4')
        val = sorted([group['name'] for group in grouping])
        ans = ['A', 'B', 'C', 'D', 'E']
        self.assertEquals(ans, val)

    def test_can_process_simple_structure(self):
        val = self.grouped_names('1G59')
        ans = [['B'], ['D']]
        self.assertEquals(ans, val)

    def test_it_can_process_a_triplet(self):
        val = self.grouped_names('4V42')
        ans = [['AA'], ['AB', 'A1'], ['AC', 'A1'], ['AD'], ['BA'], ['BB']]
        self.assertEquals(ans, val)


class SimilarSpeciesTest(StageTest):
    loader_class = Grouper

    def test_says_true_if_same(self):
        g1 = basic_data(1, species=1)
        g2 = basic_data(2, species=1)
        assert self.loader.are_similar_species(g1, g2) is True

    def test_says_different_if_disagree(self):
        g1 = basic_data(1, species=1)
        g2 = basic_data(2, species=2)
        self.assertFalse(self.loader.are_similar_species(g1, g2))

    def test_ignores_difference_if_synethic(self):
        g1 = basic_data(1, species=1)
        g2 = basic_data(2, species=32630)
        self.assertTrue(self.loader.are_similar_species(g1, g2))

    def test_ignores_difference_if_none(self):
        g1 = basic_data(1, species=1)
        g2 = basic_data(2, species=None)
        self.assertTrue(self.loader.are_similar_species(g1, g2))


class DiscrepancyTest(StageTest):
    loader_class = Grouper

    def setUp(self):
        super(DiscrepancyTest, self).setUp()
        self.disc = {
            1: {2: 0.4, 4: 1},
            2: {1: 0.4, 2: 0.0, 3: 5.0, 4: 5.0},
            3: {3: 0, 5: 0.1},
            4: {},
        }
        self.good_disc = ft.partial(self.loader.has_good_discrepancy,
                                    self.disc)

    def test_it_says_true_if_none_for_first(self):
        g1 = basic_data(2)
        g2 = basic_data(1)
        assert self.good_disc(g1, g2) is True

    def test_it_says_false_if_none_for_second(self):
        g1 = basic_data(1)
        g2 = basic_data(3)
        assert self.good_disc(g1, g2) is False

    def test_it_says_false_if_over_cutoff(self):
        g1 = basic_data(1)
        g2 = basic_data(4)
        assert self.good_disc(g1, g2) is False

    def test_it_says_true_if_under_cutoff(self):
        g1 = basic_data(1)
        g2 = basic_data(2)
        assert self.good_disc(g1, g2) is True

    def test_says_true_if_too_big_resolution(self):
        g1 = basic_data(1, resolution=5.0)
        g2 = basic_data(4)
        assert self.good_disc(g1, g2) is True
        assert self.good_disc(g2, g1) is True
        assert self.good_disc(g1, g1) is True

    def test_says_true_if_no_resolution(self):
        g1 = basic_data(1, resolution=None)
        g2 = basic_data(4, resolution=1.0)
        assert self.good_disc(g1, g2) is True
        assert self.good_disc(g2, g1) is True
        assert self.good_disc(g1, g1) is True

    def test_says_true_if_too_short(self):
        g1 = basic_data(1, length=2)
        g2 = basic_data(4, length=10)
        assert self.good_disc(g1, g2) is True
        assert self.good_disc(g2, g1) is True
        assert self.good_disc(g1, g1) is True


class AlignmentTest(StageTest):
    loader_class = Grouper

    def setUp(self):
        super(AlignmentTest, self).setUp()
        self.align = {
            1: {2: True, 3: False},
            2: {1: True, 3: False},
            3: {5: False, 1: False}
        }
        self.method = ft.partial(self.loader.has_good_alignment, self.align)

    def test_is_false_if_bad_alignment(self):
        g1 = {'id': 1, 'db_id': 1, 'chains': [{'db_id': 1}]}
        g2 = {'id': 1, 'db_id': 3, 'chains': [{'db_id': 3}]}
        self.assertFalse(self.method(g1, g2))

    def test_is_true_if_good_alignment(self):
        g1 = {'id': 1, 'db_id': 1, 'chains': [{'db_id': 1}]}
        g2 = {'id': 1, 'db_id': 2, 'chains': [{'db_id': 2}]}
        self.assertTrue(self.method(g1, g2))


class AreEquivlantTests(StageTest):
    loader_class = Grouper

    @pytest.mark.skip()
    def test_are_different_are_different_species(self):
        pass

    @pytest.mark.skip()
    def test_are_different_if_bad_align(self):
        pass

    @pytest.mark.skip()
    def test_are_different_if_bad_disc(self):
        pass

    def test_always_joins_if_given_hardcoded_chains(self):
        g1 = {'id': '1S72|1|0'}
        g2 = {'id': '1FG0|1|A'}
        g3 = {'id': '1FFZ|1|A'}
        self.assertTrue(self.loader.are_equivalent({}, {}, g1, g2))
        self.assertTrue(self.loader.are_equivalent({}, {}, g1, g3))


class LoadingIfeTest(StageTest):
    loader_class = Grouper

    def setUp(self):
        super(LoadingIfeTest, self).setUp()
        self.data = self.loader.ifes('4V9Q')

    def test_loads_all_ife_in_pdb(self):
        self.assertEquals(10, len(self.data))

    def test_assigns_species_by_first_chain(self):
        assert self.data[3]['id'] == '4V9Q|1|BV'
        self.assertEquals(512, self.data[3]['species'])

    def test_sets_length_to_first(self):
        self.assertEquals(77, self.data[3]['length'])

    def test_sets_resolution_to_first(self):
        self.assertEquals(3.4, self.data[3]['resolution'])

    def test_sets_bps(self):
        self.assertEquals(1245, self.data[0]['bp'])

    def test_sets_db_id_to_first_chain_db_id(self):
        val = self.data[0]
        print(val)
        self.assertEquals(val['db_id'], val['chains'][0]['db_id'])
        self.assertEquals(1527, val['db_id'])

    def test_sorts_ife_chains_by_stored_index(self):
        ife = self.data[3]
        self.assertEquals(ife['id'], '4V9Q|1|BV')
        val = [c['name'] for c in ife['chains']]
        ans = ['BV', 'BX']
        self.assertEquals(ans, val)

    def test_builds_correct_ife_to_chain_mappings(self):
        val = {}
        for ife in self.data:
            val[ife['id']] = []
            for chain in ife['chains']:
                val[ife['id']].append(chain['name'])
        ans = {
            '4V9Q|1|AA': ['AA'],
            '4V9Q|1|AB': ['AB'],
            '4V9Q|1|BA': ['BA'],
            '4V9Q|1|BV': ['BV', 'BX'],
            '4V9Q|1|BW': ['BW'],
            '4V9Q|1|CA': ['CA'],
            '4V9Q|1|CB': ['CB'],
            '4V9Q|1|DA': ['DA'],
            '4V9Q|1|DV': ['DV', 'DX'],
            '4V9Q|1|DW': ['DW'],
        }
        self.assertEquals(ans, val)

    def test_does_not_duplicate_ifes(self):
        ifes = self.loader.ifes('1A34')
        val = set(ife['db_id'] for ife in ifes)
        self.assertEquals(1, len(val))

    @pytest.mark.skip(reason="No data yet")
    def test_it_only_loads_ife_with_normalized_seq(self):
        pass


class PairsTest(StageTest):
    loader_class = Grouper

    def setUp(self):
        super(PairsTest, self).setUp()
        self.disc = {
            0: {1: 0.4, 3: 1, 4: 0.2},
            2: {2: 0, 4: 0.1}
        }
        self.align = {
            0: {0: True, 1: True, 2: False, 3: False, 4: True},
            1: {},
            2: {0: False, 1: False, 2: True, 3: True, 4: False},
            3: {},
            4: {},
        }
        self.chains = [
            basic_data(0, species=None),
            basic_data(1, species=1),
            basic_data(2, species=32630),
            basic_data(3, species=1),
            basic_data(4, species=1),
        ]

    def test_will_create_all_valid_pairs(self):
        assert list(self.loader.pairs(self.chains, self.align, self.disc)) == [
            (self.chains[0], self.chains[1]),
            (self.chains[0], self.chains[4]),
        ]

    def test_will_create_pairs_despite_bad_second(self):
        self.align[1][-1] = False
        assert list(self.loader.pairs(self.chains, self.align, self.disc)) == [
            (self.chains[0], self.chains[1]),
            (self.chains[0], self.chains[4]),
        ]


class ValidIfeTest(StageTest):
    loader_class = Grouper

    def test_it_will_reject_0_length_ife(self):
        self.assertFalse(self.loader.valid_ife({'id': '1', 'length': 0}))

    def test_it_will_accept_1_length_ife(self):
        self.assertTrue(self.loader.valid_ife({'id': '1', 'length': 1}))


class SpeciesSplittingTest(StageTest):
    loader_class = Grouper

    def group(self, *pdbs):
        groups = []
        self.loader.use_discrepancy = False
        for group in self.loader(pdbs):
            groups.append(sorted(set(ife['id'] for ife in group['members'])))
        self.loader.use_discrepancy = True
        return sorted(groups)

    def test_it_leaves_a_group_with_small_members_alone(self):
        assert self.group('4A3J', '4A3G', '3J9M') == [
            ['3J9M|1|A'],
            ['3J9M|1|AA'],
            ['3J9M|1|B'],
            sorted(['4A3J|1|P', '4A3G|1|P', '3J9M|1|u'])
        ]

    def test_it_splits_large_by_species(self):
        assert self.group('4V9Q', '3CW5', '4TUE') == [
            ['3CW5|1|A'],
            sorted(['4TUE|1|QA', '4TUE|1|XA', '4V9Q|1|BA', '4V9Q|1|DA']),  # Tt SSU
            sorted(['4TUE|1|QV', '4TUE|1|XV']),  # Limbella pachyloma tRNA
            ['4TUE|1|RA', '4TUE|1|YA', '4V9Q|1|AA', '4V9Q|1|CA'],
            ['4TUE|1|RB', '4TUE|1|YB', '4V9Q|1|AB', '4V9Q|1|CB'],
            ['4V9Q|1|BV', '4V9Q|1|BW', '4V9Q|1|DV', '4V9Q|1|DW'],
        ]

    @pytest.mark.skip()
    def test_it_puts_synthenic_with_largest_group(self):
        pass


class GroupingTest(StageTest):
    loader_class = Grouper

    def group(self, *pdbs):
        groups = []
        for group in self.loader(pdbs):
            groups.append(set(ife['id'] for ife in group['members']))
        return groups

    @pytest.mark.skip()
    def test_it_will_complain_if_no_ifes_found(self):
        pass

    @pytest.mark.skip()
    def test_will_filter_out_invalid_ifes(self):
        pass

    @pytest.mark.skip()
    def test_will_create_groups(self):
        pass

    @pytest.mark.skip()
    def test_will_assign_rank_to_all_grouped_members(self):
        pass

    def test_it_will_enforce_splitting_by_species(self):
        groups = self.group('4V9Q', '3CW5', '4TUE')
        assert set(['4TUE|1|QV', '4TUE|1|XV']) in groups
        assert set(['4TUE|1|QA', '4TUE|1|XA', '4V9Q|1|BA', '4V9Q|1|DA']) in groups
