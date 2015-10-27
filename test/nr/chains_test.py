from test import StageTest

from pymotifs import core
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
        self.mapping = {
            'NR_1.5_01181.1': 1,
            'NR_1.5_08345.1': 510,
            'NR_1.5_23793.1': 1776
        }
        self.grouping = [
            {
                'name': {'full': 'NR_1.5_01181.1'},
                'release': '0.1',
                'members': [{'id': 1, 'rank': 3}]
            },
            {
                'name': {'full': 'NR_1.5_08345.1'},
                'release': '0.1',
                'members': [{'id': 2, 'rank': 2}]
            },
            {
                'name': {'full': 'NR_1.5_23793.1'},
                'release': '0.1',
                'members': [{'id': 3, 'rank': 0}]
            }
        ]

    def test_maps_all_entries_in_grouping(self):
        val = self.loader.chains(self.grouping, self.mapping)
        self.assertEquals(len(self.grouping), len(val))

    def test_fails_without_a_mapping(self):
        self.assertRaises(core.InvalidState, self.loader.chains,
                          self.grouping, {})

    def test_fails_without_a_grouping(self):
        self.assertRaises(core.InvalidState, self.loader.chains,
                          {}, self.mapping)

    def test_fails_if_entry_is_not_in_mapping(self):
        copy = dict(self.grouping[2])
        copy['name']['full'] = 'bob'
        copy['release'] = 1.0
        self.grouping.append(copy)
        self.assertRaises(core.InvalidState, self.loader.chains,
                          self.grouping, self.mapping)

    def test_assigns_correct_data(self):
        val = sorted(self.loader.chains(self.grouping, self.mapping))
        ans = sorted([
            {
                'ife_id': 1,
                'nr_class_id': 1,
                'nr_release_id': '0.1',
                'rank': 3,
                'rep': False
            },
            {
                'ife_id': 2,
                'nr_class_id': 510,
                'nr_release_id': '0.1',
                'rank': 2,
                'rep': False
            },
            {
                'ife_id': 3,
                'nr_class_id': 1776,
                'nr_release_id': '0.1',
                'rank': 0,
                'rep': True
            },
        ])
        self.assertEquals(ans, val)
