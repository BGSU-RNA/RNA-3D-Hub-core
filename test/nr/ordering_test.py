import pytest

from pymotifs import core
from pymotifs.nr.ordering import Loader

from test import StageTest


class ProcessTest(StageTest):
    loader_class = Loader

    def test_it_can_get_all_class_ids(self):
        val = self.loader.to_process([])
        assert len(val) == 1090


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_no_data(self):
        assert self.loader.has_data(-1) is False

    @pytest.mark.skip()
    def test_knows_if_has_data(self):
        assert self.loader.has_data(1) is True


class MembersTest(StageTest):
    loader_class = Loader

    def test_can_load_all_members_for_a_class(self):
        assert self.loader.members(9) == [
            ('1UTD|1|0', 12),
            ('1UTD|1|2', 13),
            ('1UTD|1|3', 14),
            ('1UTD|1|4', 15),
            ('1UTD|1|5', 16),
            ('1UTD|1|6', 17),
            ('1UTD|1|7', 18),
            ('1UTD|1|8', 19),
            ('1UTD|1|9', 20),
            ('1UTD|1|Z', 21),
            ('1UTD|1|1', 22),
        ]

    def test_it_will_fail_if_no_members(self):
        with pytest.raises(core.InvalidState):
            self.loader.members(-1)


class DistancesTest(StageTest):
    loader_class = Loader

    def test_it_can_load_all_distances(self):
        members = self.loader.members(9)
        distances = self.loader.distances(9, members)
        assert len(distances) == len(members)
        assert distances['1UTD|1|0'] == {
            '1UTD|1|1': 1.11285,
            '1UTD|1|2': 0.130253,
            '1UTD|1|3': 0.0671009,
            '1UTD|1|4': 0.0800629,
            '1UTD|1|5': 0.0636083,
            '1UTD|1|6': 0.0980047,
            '1UTD|1|7': 0.0930357,
            '1UTD|1|8': 0.0721325,
            '1UTD|1|9': 0.136114,
            '1UTD|1|Z': 0.079315,
        }

    def test_it_will_raise_skip_if_no_members(self):
        with pytest.raises(core.Skip):
            self.loader.distances(-1, [])
