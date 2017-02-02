import pytest

from test import StageTest

from pymotifs.nr.groups.simplified import Grouper

from pymotifs.nr.representatives import Naive
from pymotifs.nr.representatives import Increase
from pymotifs.nr.representatives import bp_per_nt
from pymotifs.nr.representatives import ParentIncrease
from pymotifs.nr.representatives import QualityMetrics


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
