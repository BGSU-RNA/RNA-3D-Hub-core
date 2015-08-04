from test import StageTest

from pymotifs.nr.classes import Loader


class LoadReleaseTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(LoadReleaseTest, self).setUp()
        self.data = self.loader.load_release('0.1', cutoff='4.0')

    def test_loads_all(self):
        self.assertEqual(339, len(self.data))

    def test_loads_correct_data(self):
        group = filter(lambda g: g['name']['handle'] == '00396', self.data)
        self.assertEqual(1, len(group))
        group = group[0]
        name = {'handle': '00396', 'version': 1, 'full_name': 'NR_4.0_00396.1'}
        self.assertTrue('class_id' in group['name'])
        group['name'].pop('class_id')
        self.assertEqual(name, group['name'])
        self.assertEqual(5, len(group['members']))

    def test_loads_correct_members(self):
        group = filter(lambda g: g['name']['handle'] == '00396', self.data)
        self.assertEqual(1, len(group))
        group = group[0]
        ans = sorted([{'id': '2HOJ|A'}, {'id': '2HOK|A'}, {'id': '2HOL|A'},
                      {'id': '2HOM|A'}, {'id': '2HOO|A'}])
        self.assertEqual(ans, sorted(group['members']))

    def test_it_loads_ignoring_structure_resolution(self):
        group = filter(lambda g: g['name']['handle'] == '00834', self.data)
        self.assertEqual(1, len(group))
        self.assertEqual(49, len(group[0]['members']))


class WithinCutoffTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(WithinCutoffTest, self).setUp()
        self.chains = [
            {'resolution': 1.0, 'rank': 3},
            {'resolution': 2.0, 'rank': 0},
            {'resolution': 3.0, 'rank': 1},
            {'resolution': 3.5, 'rank': 4},
            {'resolution': 5.0, 'rank': 5},
            {'resolution': 4.0, 'rank': 6}
        ]

    def test_it_can_select_all_with_all_cutoff(self):
        val = self.loader.within_cutoff(self.chains, 'all')
        ans = sorted(self.chains, key=lambda c: c['rank'])
        self.assertEqual(ans, val)

    def test_it_can_select_all_below_cutoff(self):
        val = self.loader.within_cutoff(self.chains, '4.0')
        ans = [self.chains[1], self.chains[2], self.chains[0], self.chains[3],
               self.chains[5]]
        ans[0]['rank'] = 0
        ans[1]['rank'] = 1
        ans[2]['rank'] = 2
        ans[3]['rank'] = 3
        ans[4]['rank'] = 4
        self.assertEqual(ans, val)
