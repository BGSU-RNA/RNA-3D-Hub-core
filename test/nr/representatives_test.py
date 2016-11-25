import pytest

from test import StageTest

from pymotifs.nr.groups.simplified import Grouper

from pymotifs.nr.representatives import Naive
from pymotifs.nr.representatives import Increase
from pymotifs.nr.representatives import bp_per_nt
from pymotifs.nr.representatives import ParentIncrease


class SortingChainsTest(StageTest):

    def test_builds_key_with_nt_per_bp_and_id(self):
        chain1 = {'bp': 10, 'length': 20, 'id': 'a'}
        chain2 = {'bp': 20, 'length': 20, 'id': 'b'}
        self.assertEquals(bp_per_nt(chain1), (0.5, None, 'a'))
        self.assertEquals(bp_per_nt(chain2), (1, None, 'b'))

    def test_builds_key_with_0_if_one_is_0(self):
        chain1 = {'bp': 0, 'length': 10, 'id': 'c'}
        chain2 = {'bp': 10, 'length': 0, 'id': 'd'}
        chain3 = {'bp': 0, 'length': 0, 'id': 'e'}
        self.assertEquals(bp_per_nt(chain1), (0, None, 'c'))
        self.assertEquals(bp_per_nt(chain2), (0, None, 'd'))
        self.assertEquals(bp_per_nt(chain3), (0, None, 'e'))

    def test_uses_inverse_resolution(self):
        chain1 = {'bp': 0, 'length': 10, 'id': 'c', 'resolution': 10}
        self.assertEquals(bp_per_nt(chain1), (0, -10, 'c'))


class NaiveBestTest(StageTest):
    loader_class = Naive

    def best(self, *chains):
        return self.loader({'members': chains, 'parent': []})

    def test_it_can_deal_with_0_bp_or_length_chains(self):
        val = self.best({'bp': 13, 'length': 10, 'id': 'c'},
                        {'bp': 0, 'length': 10, 'id': 'a'},
                        {'bp': 50, 'length': 0, 'id': 'b'})
        ans = {'bp': 13, 'length': 10, 'id': 'c'}
        self.assertEquals(ans, val)

    def test_gets_by_bp_per_nt(self):
        val = self.best({'bp': 13, 'length': 10, 'id': 'c'},
                        {'bp': 10, 'length': 10, 'id': 'a'},
                        {'bp': 50, 'length': 10, 'id': 'b'})
        ans = {'bp': 50, 'length': 10, 'id': 'b'}
        self.assertEquals(ans, val)

    def test_it_tiebreaks_on_pdb(self):
        val = self.best({'bp': 10, 'length': 10, 'id': 'c'},
                        {'bp': 10, 'length': 10, 'id': 'a'},
                        {'bp': 10, 'length': 10, 'id': 'b'})
        ans = {'bp': 10, 'length': 10, 'id': 'c'}
        self.assertEquals(ans, val)

    def test_if_all_have_0_bp_uses_resolution(self):
        val = self.best({'bp': 0, 'length': 10, 'id': 'c', 'resolution': 4.0},
                        {'bp': 0, 'length': 10, 'id': 'a', 'resolution': 3.0},
                        {'bp': 0, 'length': 10, 'id': 'b', 'resolution': 2.0})
        ans = {'bp': 0, 'length': 10, 'id': 'b', 'resolution': 2.0}
        self.assertEquals(ans, val)

    @pytest.mark.skip()
    def test_it_ignores_any_set_parent(self):
        val = self.loader({
            'parent': [{'representative': {'bp': 100, 'length': 10}}],
            'members': [{'bp': 13, 'length': 10, 'id': 'c'},
                        {'bp': 0, 'length': 10, 'id': 'a'},
                        {'bp': 50, 'length': 0, 'id': 'b'}],
        })
        ans = {'bp': 13, 'length': 10, 'id': 'c'}
        assert ans == val


class IncreaseCandidatesTest(StageTest):
    loader_class = Increase

    def candidates(self, *chains):
        default = {'bp': 10, 'length': 20, 'id': 'a'}
        return self.loader.candidates(default, chains)

    def test_it_can_get_all_longer_ifes(self):
        val = self.candidates({'bp': 10, 'length': 10, 'id': 'b'},
                              {'bp': 9, 'length': 100, 'id': 'e'},
                              {'bp': 10, 'length': 30, 'id': 'c'},
                              {'bp': 10, 'length': 20, 'id': 'd'})
        ans = [{'bp': 10, 'length': 30, 'id': 'c'}]
        self.assertEquals(ans, val)

    def test_it_can_get_all_with_more_bp(self):
        val = self.candidates({'bp': 10, 'length': 20, 'id': 'b'},
                              {'bp': 100, 'length': 1, 'id': 'e'},
                              {'bp': 13, 'length': 20, 'id': 'c'},
                              {'bp': 0, 'length': 20, 'id': 'd'})
        ans = [{'bp': 13, 'length': 20, 'id': 'c'}]
        self.assertEquals(ans, val)

    def test_it_can_sort_results_correctly(self):
        val = self.candidates({'bp': 11, 'length': 20, 'id': 'b'},
                              {'bp': 0, 'length': 20, 'id': 'd'},
                              {'bp': 13, 'length': 20, 'id': 'c'})
        ans = [{'bp': 13, 'length': 20, 'id': 'c'},
               {'bp': 11, 'length': 20, 'id': 'b'}]
        self.assertEquals(ans, val)


class ComputingIncreaseTest(StageTest):
    loader_class = Increase

    def increase(self, chain1, chain2, key='a'):
        return self.loader.increase(chain1, chain2, key)

    def test_it_computes_percent_increase(self):
        val = self.increase({'a': 10}, {'a': 1})
        self.assertEquals(val, 900.0)

    def test_it_can_handle_stays_same(self):
        val = self.increase({'a': 10}, {'a': 10})
        self.assertEquals(val, 0.0)

    def test_it_can_handle_decrease(self):
        val = self.increase({'a': 5}, {'a': 10})
        self.assertEquals(val, -50.0)

    def test_it_gives_100_for_given_0(self):
        self.assertEquals(self.increase({'a': 0}, {'a': 10}), -100)
        self.assertEquals(self.increase({'a': 10}, {'a': 0}), 100)

    def test_given_both_0_it_gives_0(self):
        self.assertEquals(self.increase({'a': 0}, {'a': 0}), 0)


class FilteringGroupTest(StageTest):
    loader_class = Increase

    def test_it_can_filter_a_group_to_correct_method(self):
        val = self.loader.filter_group({
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 3, 'method': 'other'}],
            'parent': []
        })
        assert val == {
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'}],
            'parent': []
        }

    def test_it_can_filter_a_group_to_given_methods(self):
        val = self.loader.filter_group({
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 3, 'method': 'other1'},
                        {'id': 4, 'method': 'other2'}],
            'parent': []
        }, methods=set(['other1', 'other2']))
        assert val == {
            'members': [{'id': 3, 'method': 'other1'},
                        {'id': 4, 'method': 'other2'}],
            'parent': []
        }

    def test_it_will_keep_all_if_none_match_method(self):
        val = self.loader.filter_group({
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 3, 'method': 'other'}],
            'parent': []
        }, methods=set(['other1', 'other2']))
        assert val == {
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 3, 'method': 'other'}],
            'parent': []
        }

    def test_it_will_use_given_parent(self):
        val = self.loader.filter_group({
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 3, 'method': 'other'}],
            'parent': [{'id': 3}],
        })
        assert val == {
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'}],
            'parent': [{'id': 3}],
        }

    def test_it_will_set_parent_to_empty_if_none(self):
        val = self.loader.filter_group({
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 3, 'method': 'other'}],
        })
        assert val == {
            'members': [{'id': 1, 'method': 'X-RAY DIFFRACTION'},
                        {'id': 2, 'method': 'X-RAY DIFFRACTION'}],
            'parent': []
        }



class CyroEmTest(StageTest):
    loader_class = Increase

    def rep(self, group):
        return self.loader({'members': group})

    def test_it_prefers_xray_to_cyro(self):
        chains = [
            {'bp': 10, 'length': 10, 'resolution': 4, 'id': 'b', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 10, 'length': 10, 'resolution': 4, 'id': 'a', 'method': 'X-RAY DIFFRACTION'},
        ]
        assert self.rep(chains) == chains[1]

    def test_it_will_use_xray_with_worse_resolution_than_cyro(self):
        chains = [
            {'bp': 10, 'length': 10, 'resolution': 2, 'id': 'b', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 10, 'length': 10, 'resolution': 3, 'id': 'a', 'method': 'X-RAY DIFFRACTION'},
            {'bp': 20, 'length': 10, 'resolution': 3, 'id': 'c', 'method': 'X-RAY DIFFRACTION'},
        ]
        assert self.rep(chains) == chains[2]

    def test_it_will_use_xray_over_better_bp_nt_cyro(self):
        chains = [
            {'bp': 40, 'length': 20, 'resolution': 2, 'id': 'b', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 10, 'length': 10, 'resolution': 3, 'id': 'a', 'method': 'X-RAY DIFFRACTION'},
            {'bp': 10, 'length': 10, 'resolution': 4, 'id': 'c', 'method': 'X-RAY DIFFRACTION'},
        ]
        assert self.rep(chains) == chains[1]

    def test_it_will_use_cryo_if_only_cyro(self):
        chains = [
            {'bp': 10, 'length': 10, 'resolution': 3, 'id': 'b', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 20, 'length': 20, 'resolution': 3, 'id': 'a', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 20, 'length': 12, 'resolution': 3, 'id': 'c', 'method': 'ELECTRON MICROSCOPY'},
        ]
        assert self.rep(chains) == chains[2]


class PickingRepresentativeTest(StageTest):
    loader_class = Increase

    def group(self, *args, **kwargs):
        grouper = Grouper(self.loader.config, self.loader.session)
        members = []
        for pdb, chains in args:
            ifes = grouper.ifes(pdb)
            members.extend(ife for ife in ifes if ife['name'] in chains)
        return members

    def rep(self, *args, **kwargs):
        group = self.group(*args, **kwargs)
        return self.loader({'members': group})

    def test_4MGM_4MGN(self):
        val = self.rep(('4MGM', ('A', 'B')), ('4MGN', ('B', 'D')),
                       id='NR_4.0_13428.1')
        assert val['id'] == '4MGN|1|B'

    def test_1GID(self):
        val = self.rep(('1GID', ('A', 'B')), ('1X8W', ('A', 'B', 'C', 'D')),
                       ('1GRZ', ('A', 'B')), id='NR_4.0_86492.1')
        assert val['id'] == '1X8W|1|A'

    @pytest.mark.xfail(reason="No cutoffs for small chains yet."
                       " Selects: 4OAU|1|A")
    def test_4OAU(self):
        val = self.rep(('4OAU', 'A'), ('3S4G', 'B'), ('3GPQ', ('E', 'F')),
                       id='NR_all_14757.1')
        assert val['id'] == '3GPQ|1|E'


class ParentIncreaseTest(StageTest):
    loader_class = ParentIncrease

    def test_it_will_select_parent_representative_if_exists(self):
        val = self.loader.initial_representative({
            'members': [{'id': 'A', 'bp': 10, 'length': 10},
                        {'id': 'B', 'bp': 8, 'length': 10}],
            'parent': [{'representative': {'id': 'B', 'bp': 8, 'length': 10}}]
        })
        assert val == {'id': 'B', 'bp': 8, 'length': 10}

    def test_it_will_use_naive_if_no_parent(self):
        val = self.loader.initial_representative({
            'members': [{'id': 'A', 'bp': 10, 'length': 10},
                        {'id': 'B', 'bp': 8, 'length': 10}],
            'parent': []
        })
        assert val == {'id': 'A', 'bp': 10, 'length': 10}

    def test_will_not_use_parent_if_not_current_member(self):
        val = self.loader.initial_representative({
            'members': [{'id': 'A', 'bp': 10, 'length': 10},
                        {'id': 'C', 'bp': 8, 'length': 10}],
            'parent': [{'representative': {'id': 'B', 'bp': 8, 'length': 10}}]
        })
        assert val == {'id': 'A', 'bp': 10, 'length': 10}

    def test_if_too_many_parents_uses_naive(self):
        val = self.loader.initial_representative({
            'members': [{'id': 'A', 'bp': 10, 'length': 10},
                        {'id': 'C', 'bp': 8, 'length': 10}],
            'parent': [
                {'representative': {'id': 'D', 'bp': 7, 'length': 9}},
                {'representative': {'id': 'B', 'bp': 8, 'length': 10}}
            ]
        })
        assert val == {'id': 'A', 'bp': 10, 'length': 10}

    @pytest.mark.skip()
    def test_it_requires_percent_increase_over_parent(self):
        pass

    @pytest.mark.skip()
    def test_it_will_not_change_with_small_increase(self):
        pass
