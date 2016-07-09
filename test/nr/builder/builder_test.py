import pytest

from test import StageTest

from pymotifs.core import InvalidState
from pymotifs.constants import RESOLUTION_GROUPS as GROUPS
from pymotifs.nr.builder import Builder
from pymotifs.nr.builder import Known


class GroupingTest(StageTest):
    loader_class = Builder

    def group(self, pdbs):
        grouped = self.loader.group(pdbs)
        for group in grouped:
            for member in group['members']:
                del member['db_id']
                for chain in member['chains']:
                    del chain['db_id']
        return grouped

    def test_complains_if_nothing_to_group(self):
        with pytest.raises(InvalidState):
            self.loader.group([])

    def test_will_group_given_pdbs(self):
        val = self.group(['1GID'])
        assert val == [{
            'members': [
                {'length': 158, 'bp': 76, 'name': 'B', 'species': 32630, 'resolution': 2.5, 'chains': [{u'is_accompanying': 0, u'is_integral': 1, 'name': 'B', u'sequence': 'GAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUC', 'pdb': '1GID', u'length': 158L, 'bp': 76L, 'species': 32630L, u'resolution': 2.5, 'id': '1GID|1|B'}], 'id': '1GID|1|B', 'rank': 0, 'pdb': '1GID'},
                {'length': 158, 'bp': 73, 'name': 'A', 'species': 32630, 'resolution': 2.5, 'chains': [{u'is_accompanying': 0, u'is_integral': 1, 'name': 'A', u'sequence': 'GAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUC', 'pdb': '1GID', u'length': 158L, 'bp': 73L, 'species': 32630L, u'resolution': 2.5, 'id': '1GID|1|A'}], 'id': '1GID|1|A', 'rank': 1, 'pdb': '1GID'}
            ]
        }]


class FilteringTest(StageTest):
    loader_class = Builder

    def setUp(self):
        super(FilteringTest, self).setUp()
        self.data = {
            'id': 1,
            'members': [
                {'id': 'a', 'resolution': 3},
                {'id': 'b', 'resolution': 1.0},
                {'id': 'c', 'resolution': 3},
                {'id': 'd', 'resolution': 2.0},
                {'id': 'e', 'resolution': 10.0},
                {'id': 'f', 'resolution': None},
            ]
        }

    def test_can_filter_all_within_cutoff(self):
        filtered = self.loader.within_cutoff(self.data, 2.0)
        assert filtered == {
            'id': 1,
            'members': [
                {'id': 'b', 'resolution': 1.0, 'rank': 0},
                {'id': 'd', 'resolution': 2.0, 'rank': 1},
            ]
        }

    def test_will_return_all_given_all(self):
        assert self.loader.within_cutoff(self.data, 'all') == self.data

    def test_will_use_float_cutoff(self):
        filtered = self.loader.within_cutoff(self.data, '4.0')
        assert filtered == {
            'id': 1,
            'members': [
                {'id': 'a', 'resolution': 3, 'rank': 0},
                {'id': 'b', 'resolution': 1.0, 'rank': 1},
                {'id': 'c', 'resolution': 3, 'rank': 2},
                {'id': 'd', 'resolution': 2.0, 'rank': 3},
            ]
        }

    def test_will_return_empty_given_nothing_that_passes(self):
        filtered = self.loader.within_cutoff(self.data, '0.5')
        assert filtered == {}


class NamingTest(StageTest):
    loader_class = Builder

    def named(self, release1, release2):
        known = Known(self.loader.config, self.loader.session)
        parents = known.classes(release1, '4.0')
        current = known.classes(release2, '4.0')
        data = {}
        for group in self.loader.named(current, parents):
            data[(group['name']['handle'], group['name']['version'])] = group
        return data

    def test_it_will_name_given_groups(self):
        val = self.named('1.0', '2.0')
        assert val
        for (handle, version) in val.keys():
            assert handle
            assert version


class ParentCountsTest(StageTest):
    loader_class = Builder

    def counts(self, release1, release2):
        known = Known(self.loader.config, self.loader.session)
        parents = known.classes(release1, 'all')
        groups = known.classes(release2, 'all')
        return self.loader.cutoff_counts(parents, groups)

    def test_can_compute_correct_parent_counts_for_one_cutoff(self):
        assert self.counts('1.0', '2.0') == {
            'ifes': {
                'added': 39,
                'removed': 3,
                'unchanged': 117,
            },
            'classes': {
                'added': 14,
                'removed': 4,
                'updated': 4,
                'unchanged': 32,
            },
            'pdbs': {
                'added': 9,
                'removed': 2,
            }
        }


class BuildingTest(StageTest):
    loader_class = Builder

    def test_can_build_reasonable_sets(self):
        val = self.loader(["124D", "157D"], '1.0', '2.0')
        assert val
        assert sorted(val.keys()) == ['groups', 'parent', 'parent_counts',
                                      'release']
        assert val['release'] == '2.0'
        assert val['parent'] == '1.0'

        parent_counts = sorted(val['parent_counts'],
                               key=lambda c: GROUPS.index(c['cutoff']))
        assert parent_counts == [
            {'classes': {'added': 0, 'removed': 3, 'unchanged': 1, 'updated': 0},
             'ifes': {'added': 0, 'removed': 4, 'unchanged': 1},
             'pdbs': {'added': 0, 'removed': 3},
             'cutoff': '2.0'},
            {'classes': {'added': 0, 'removed': 8, 'unchanged': 1, 'updated': 0},
             'ifes': {'added': 0, 'removed': 25, 'unchanged': 1},
             'pdbs': {'added': 0, 'removed': 8},
             'cutoff': '2.5'},
            {'classes': {'added': 0, 'removed': 17, 'unchanged': 1, 'updated': 0},
             'ifes': {'added': 0, 'removed': 65, 'unchanged': 1},
             'pdbs': {'added': 0, 'removed': 17},
             'cutoff': '3.0'},
            {'classes': {'added': 0, 'removed': 28, 'unchanged': 1, 'updated': 0},
             'ifes': {'added': 0, 'removed': 97, 'unchanged': 1},
             'pdbs': {'added': 0, 'removed': 30},
             'cutoff': '3.5'},
            {'classes': {'added': 0, 'removed': 31, 'unchanged': 1, 'updated': 0},
             'ifes': {'added': 0, 'removed': 108, 'unchanged': 1},
             'pdbs': {'added': 0, 'removed': 33},
             'cutoff': '4.0'},
            {'classes': {'added': 0, 'removed': 36, 'unchanged': 1, 'updated': 0},
             'ifes': {'added': 0, 'removed': 113, 'unchanged': 1},
             'pdbs': {'added': 0, 'removed': 34},
             'cutoff': '20.0'},
            {'classes': {'added': 1, 'removed': 39, 'unchanged': 1, 'updated': 0},
             'ifes': {'added': 1, 'removed': 119, 'unchanged': 1},
             'pdbs': {'added': 1, 'removed': 37},
             'cutoff': 'all'},
        ]

        assert len(val['groups']) == 8
