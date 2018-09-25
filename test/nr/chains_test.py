from test import StageTest

from pymotifs.nr.chains import Loader


# class LoadingIdsTest(StageTest):
#     loader_class = Loader

#     def test_gets_all_ids(self):
#         names = ['NR_1.5_01181.1', 'NR_1.5_08345.1', 'NR_1.5_23793.1',
#                  'NR_1.5_26877.1', 'NR_1.5_31163.1']
#         val = self.loader.ids(names, '0.1')
#         self.assertEqual(len(names), len(val))

#     def test_gets_correct_ids(self):
#         names = ['NR_1.5_01181.1', 'NR_1.5_08345.1', 'NR_1.5_23793.1']
#         ans = {
#             'NR_1.5_01181.1': 1,
#             'NR_1.5_08345.1': 510,
#             'NR_1.5_23793.1': 1776
#         }
#         val = self.loader.ids(names, '0.1')
#         self.assertEqual(ans, val)

#     def test_raises_exception_if_nothing_loaded(self):
#         self.assertRaises(core.InvalidState, self.loader.ids, [], '1.0')

#     def test_raises_if_finds_too_few(self):
#         names = ['NR_1.5_01181.1', 'NR_1.5_08345.1', 'NR_1.5_23793.1',
#                  'NR_1.5_26877.1', 'bob', 'NR_1.5_31163.1']
#         self.assertRaises(core.InvalidState, self.loader.ids, names, '0.1')


class ChainsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ChainsTest, self).setUp()
        self.grouping = [
            {
                'name': {'full': 'NR_1.5_01181.1', 'class_id': 1},
                'members': [{'id': 1, 'rank': 3}, {'id': -1, 'rank': 0}],
                'representative': {'id': -1},
            },
            {
                'name': {'full': 'NR_1.5_08345.1', 'class_id': 510},
                'members': [{'id': 2, 'rank': 2}],
                'representative': {'id': 4},
            },
            {
                'name': {'full': 'NR_1.5_23793.1', 'class_id': 1776},
                'members': [{'id': 3, 'rank': 0}],
                'representative': {'id': 3}
            }
        ]

    def test_maps_all_entries_in_grouping(self):
        val = self.loader.chains('0.1', self.grouping)
        assert len(val) == 4

    def test_assigns_correct_data(self):
        assert self.loader.chains('0.1', self.grouping) == [
            {'ife_id': 1, 'nr_class_id': 1, 'nr_release_id': '0.1', 'rank': 3, 'rep': False},
            {'ife_id': -1, 'nr_class_id': 1, 'nr_release_id': '0.1', 'rank': 0, 'rep': True},
            {'ife_id': 2, 'nr_class_id': 510, 'nr_release_id': '0.1', 'rank': 2, 'rep': False},
            {'ife_id': 3, 'nr_class_id': 1776, 'nr_release_id': '0.1', 'rank': 0, 'rep': True},
        ]
