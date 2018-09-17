import pytest

from test import StageTest

from pymotifs.nr.groups.simplified import Grouper

from pymotifs.nr.representatives import Increase


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
        val = self.loader.filter_group_by_method({
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
        val = self.loader.filter_group_by_method({
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
        val = self.loader.filter_group_by_method({
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
        val = self.loader.filter_group_by_method({
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
        val = self.loader.filter_group_by_method({
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
        assert self.rep(chains) == [chains[1], chains[0]]

    def test_it_will_use_xray_with_worse_resolution_than_cyro(self):
        chains = [
            {'bp': 10, 'length': 10, 'resolution': 2, 'id': 'b', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 10, 'length': 10, 'resolution': 3, 'id': 'a', 'method': 'X-RAY DIFFRACTION'},
            {'bp': 20, 'length': 10, 'resolution': 3, 'id': 'c', 'method': 'X-RAY DIFFRACTION'},
        ]
        assert self.rep(chains) == [chains[2], chains[0], chains[1]]

    def test_it_will_use_xray_over_better_bp_nt_cyro(self):
        chains = [
            {'bp': 40, 'length': 20, 'resolution': 2, 'id': 'b', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 10, 'length': 10, 'resolution': 3, 'id': 'a', 'method': 'X-RAY DIFFRACTION'},
            {'bp': 10, 'length': 10, 'resolution': 4, 'id': 'c', 'method': 'X-RAY DIFFRACTION'},
        ]
        assert self.rep(chains) == [chains[1], chains[0], chains[2]]

    def test_it_will_use_cryo_if_only_cyro(self):
        chains = [
            {'bp': 10, 'length': 10, 'resolution': 3, 'id': 'b', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 20, 'length': 20, 'resolution': 3, 'id': 'a', 'method': 'ELECTRON MICROSCOPY'},
            {'bp': 20, 'length': 12, 'resolution': 3, 'id': 'c', 'method': 'ELECTRON MICROSCOPY'},
        ]
        assert self.rep(chains) == [chains[2], chains[0], chains[1]]


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
        return [m['id'] for m in self.loader({'members': group})]

    def test_4MGM_4MGN(self):
        val = self.rep(('4MGM', ('A', 'B')), ('4MGN', ('B', 'D')),
                       id='NR_4.0_13428.1')
        assert val == ['4MGN|1|B', '4MGM|1|A', '4MGM|1|B', '4MGN|1|D']

    def test_1GID(self):
        val = self.rep(('1GID', ('A', 'B')), ('1X8W', ('A', 'B', 'C', 'D')),
                       ('1GRZ', ('A', 'B')), id='NR_4.0_86492.1')
        assert val == [
            '1X8W|1|A',
            '1GID|1|B',
            '1GID|1|A',
            '1X8W|1|B',
            '1X8W|1|D',
            '1X8W|1|C',
            '1GRZ|1|A',
            '1GRZ|1|B',
        ]

    @pytest.mark.xfail(reason="No cutoffs for small chains yet."
                       " Selects: 4OAU|1|A")
    def test_4OAU(self):
        val = self.rep(('4OAU', 'A'), ('3S4G', 'B'), ('3GPQ', ('E', 'F')),
                       id='NR_all_14757.1')
        assert val == [
            '3GPQ|1|F',
            '3GPQ|1|E',
            '3S4G|1|B',
            '4OAU|1|A',
        ]
