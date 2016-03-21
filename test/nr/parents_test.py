import pytest

from test import StageTest

from pymotifs import core
from pymotifs.nr.parents import Loader


class ParentInfoTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ParentInfoTest, self).setUp()
        self.grouping = {
        }

    def test_it_fails_without_a_grouping(self):
        self.assertRaises(core.InvalidState, self.loader.parents, {}, {'a': 1})

    def test_it_fails_without_a_mapping(self):
        self.assertRaises(core.InvalidState, self.loader.parents, {'a': 1}, {})

    @pytest.mark.skip()
    def test_it_gets_previous_release(self):
        _, release_id = self.loader.parents(self.grouping)
        self.assertEquals('1.0', release_id)

    @pytest.mark.skip()
    def test_it_gets_all_parent_names(self):
        val, _ = self.loader.parents(self.grouping)
        ans = sorted(['NR_1.5_01181.1', 'NR_1.5_08345.1', 'NR_1.5_23793.1',
                      'NR_1.5_26877.1', 'NR_1.5_31163.1'])
        self.assertEquals(ans, sorted(val))

    @pytest.mark.skip()
    def test_fails_if_no_parent_names_found_with_previous_release(self):
        pass

    @pytest.mark.skip()
    def test_gives_no_parents_if_no_previous_release(self):
        val = self.loader.parents(self.grouping)
        ans = ([], None)
        self.assertEquals(ans, val)
