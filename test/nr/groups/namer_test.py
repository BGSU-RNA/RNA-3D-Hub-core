from test import StageTest
from nose import SkipTest

from pymotifs.nr.groups.naming import Namer


class OverlapTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(OverlapTest, self).setUp()
        self.group1 = {
            'members': [
                {'id': 'A'},
                {'id': 'B'},
                {'id': 'C'}
            ],
        }
        self.group2 = {
            'members': [
                {'id': 'A'},
                {'id': 'B'},
                {'id': 'D'},
                {'id': 'E'}
            ],
        }
        self.group3 = {
            'members': [
                {'id': 'D'},
                {'id': 'E'}
            ],
        }

    def test_it_can_determine_overlap(self):
        val = self.loader.overlap(self.group1, self.group2)
        ans = {'group': self.group2, 'intersection': set(['A', 'B'])}
        self.assertEqual(ans, val)

    def test_it_gives_empty_for_no_overlap(self):
        val = self.loader.overlap(self.group1, self.group3)
        ans = {}
        self.assertEqual(ans, val)

    def test_it_can_computes_the_min_overlap_size(self):
        overlap = self.loader.overlap(self.group1, self.group2)
        val = self.loader.overlap_size(self.group1, overlap)
        ans = 2.0 / 4.0
        self.assertEquals(ans, val)


class NewNamingTest(StageTest):
    loader_class = Namer

    def test_it_can_create_a_new_handler(self):
        naming = self.loader.new_name(1)
        self.assertEqual(5, len(naming['handle']))

    def test_sets_the_version_to_one(self):
        naming = self.loader.new_name(1)
        self.assertEqual(1, naming['version'])

    def test_sets_the_no_parent_comment(self):
        naming = self.loader.new_name(0)
        self.assertEqual('New id, no parents', naming['comment'])

    def test_sets_the_single_parent_comment(self):
        naming = self.loader.new_name(1)
        self.assertEqual('New id, 1 parent', naming['comment'])

    def test_sets_the_two_parent_comment(self):
        naming = self.loader.new_name(2)
        self.assertEqual('New id, 2 parents', naming['comment'])

    def test_sets_the_multiple_parent_comment(self):
        naming = self.loader.new_name(4)
        self.assertEqual('New id, > 2 parents', naming['comment'])


class SameNameTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(SameNameTest, self).setUp()
        self.group = {
            'handle': '00001',
            'version': 3,
            'comment': 'Something'
        }

    def test_it_copies_over_old_handle(self):
        name = self.loader.same_name(self.group)
        self.assertEqual('00001', name['handle'])

    def test_it_keeps_old_version(self):
        name = self.loader.same_name(self.group)
        self.assertEqual(3, name['version'])

    def test_it_adds_exact_comment(self):
        name = self.loader.same_name(self.group)
        self.assertEqual('Exact match', name['comment'])


class UpdatedNameTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(UpdatedNameTest, self).setUp()
        self.group = {
            'handle': '03001',
            'version': 3,
            'comment': 'Another comment'
        }

    def test_it_copies_over_old_handle(self):
        name = self.loader.updated_name(self.group, 1)
        self.assertEqual('03001', name['handle'])

    def test_it_bumps_old_version(self):
        name = self.loader.updated_name(self.group, 2)
        self.assertEqual(4, name['version'])

    def test_it_adds_comment_for_one_parent(self):
        name = self.loader.updated_name(self.group, 1)
        self.assertEqual('Updated, 1 parent', name['comment'])

    def test_it_adds_comment_for_two_parents(self):
        name = self.loader.updated_name(self.group, 2)
        self.assertEqual('Updated, 2 parents', name['comment'])


class OneParentTest(StageTest):
    loader_class = Namer

    def setUp(self):
        super(OneParentTest, self).setUp()
        self.group = {'members': [{'id': 'A'}, {'id': 'B'}, {'id': 'C'}]}
        self.exact = {
            'group': dict(self.group),
            'intersection': set(['A', 'B', 'C'])
        }
        self.exact['group']['name'] = {
            'handle': '10000',
            'comment': 'A',
            'version': 3
        }
        self.large = {
            'members': [{'id': 'A'}, {'id': 'C'}],
            'large': {'handle': '20000', 'comment': 'A', 'version': 1}
        }
        self.small = {
            'members': [{'id': 'D'}, {'id': 'C'}, {'id': 'D'}],
            'name': {'handle': '30000', 'comment': 'A', 'version': 2}
        }

    def test_it_uses_same_name_with_exact_overlap(self):
        val = self.loader.one_parent(self.group, self.exact)
        ans = {'handle': '10000', 'comment': 'Exact match', 'version': 3}
        self.assertEqual(ans, val)

    def test_it_use_updated_name_with_large_overlap(self):
        raise SkipTest()

    def test_it_gives_new_name_if_overlap_to_old_is_small(self):
        raise SkipTest()

    def test_it_gves_new_name_if_overlap_to_new_is_small(self):
        raise SkipTest()


class TwoParentTest(StageTest):
    loader_class = Namer

    def test_it_updates_if_one_old_has_large_overlap(self):
        raise SkipTest()

    def test_it_uses_new__if_both_have_small_overlap(self):
        raise SkipTest()


class ManyParentTest(StageTest):
    loader_class = Namer

    def test_it_creates_new_name(self):
        raise SkipTest()

    def test_it_ignores_overlap_to_old(self):
        raise SkipTest()
