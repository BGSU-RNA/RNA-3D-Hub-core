from test import StageTest
from nose import SkipTest

from pymotifs import core
from pymotifs.nr.parents import Loader


class ParentInfoTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ParentInfoTest, self).setUp()
        self.grouping = {
        }

    def test_it_fails_without_a_grouping(self):
        self.assertRaises(core.InvalidState, self.loader.parents, {})

    def test_it_gets_previous_release(self):
        raise SkipTest()
        # _, release_id = self.loader.parents(self.grouping)
        # self.assertEquals('1.0', release_id)

    def test_it_gets_all_parent_names(self):
        raise SkipTest()
        val, _ = self.loader.parents(self.grouping)
        ans = sorted(['NR_1.5_01181.1', 'NR_1.5_08345.1', 'NR_1.5_23793.1',
                      'NR_1.5_26877.1', 'NR_1.5_31163.1'])
        self.assertEquals(ans, sorted(val))

    def test_fails_if_no_parent_names_found_with_previous_release(self):
        raise SkipTest()
        # self.assertRaises(core.InvalidState, self.loader.parent_info

    def test_gives_no_parents_if_no_previous_release(self):
        raise SkipTest()
        val = self.loader.parents(self.grouping)
        ans = ([], None)
        self.assertEquals(ans, val)
