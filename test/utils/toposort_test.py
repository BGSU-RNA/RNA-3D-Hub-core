import pytest

from pymotifs.utils import toposort as topo

from unittest import TestCase


class LevelsTest(TestCase):
    def test_it_can_return_all_levels(self):
        data = {
            1: set([3, 4]),
            2: set([3]),
            3: set([4]),
        }
        assert len(list(topo.levels(data))) == 3

    def test_it_returns_correct_levels(self):
        data = {
            1: set([3, 4]),
            2: set([3]),
            3: set([4]),
        }
        assert list(topo.levels(data)) == [
            set([4]), set([3]), set([1, 2])
        ]

    def test_it_complains_about_circular(self):
        data = {
            1: set([3, 4]),
            2: set([3, 4]),
            3: set([4, 2]),
            4: set([2]),
        }
        with pytest.raises(ValueError):
            list(topo.levels(data))


class SortingTest(TestCase):
    def test_it_will_sort_by_given_sorter(self):
        data = {
            1: set([3, 4]),
            2: set([3]),
            3: set([4]),
        }
        assert topo.toposort(data, by=lambda n: n) == [4, 3, 1, 2]
        assert topo.toposort(data, by=lambda n: -n) == [4, 3, 2, 1]
