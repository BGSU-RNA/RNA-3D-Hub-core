from test import StageTest

from pymotifs.nr.classes import Loader


class LoadReleaseTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(LoadReleaseTest, self).setUp()
        self.data = self.loader.load_release('0.1', '4.0')

    def test_loads_all(self):
        self.assertEqual(339, len(self.data))

    def test_loads_correct_data(self):
        group = filter(lambda g: g['handle'] == '00396', self.data)
        self.assertEqual(1, len(group))
        self.assertEqual(1, group[0]['version'])
        self.assertEqual(5, len(group[0]['members']))

    def test_loads_correct_members(self):
        group = filter(lambda g: g['handle'] == '00396', self.data)
        group = group[0]
        ans = sorted([{'id': '2HOJ|A'}, {'id': '2HOK|A'}, {'id': '2HOL|A'},
                      {'id': '2HOM|A'}, {'id': '2HOO|A'}])
        self.assertEqual(ans, sorted(group['members']))
