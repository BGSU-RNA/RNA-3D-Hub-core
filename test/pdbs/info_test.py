from test import StageTest

from pymotifs.pdbs.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_knows_if_data_is_present(self):
        self.assertTrue(self.loader.has_data('2AW7'))

    def test_knows_if_data_is_missing(self):
        self.assertFalse(self.loader.has_data('bob'))

    def test_it_can_create_data(self):
        val = self.loader.data('2AW7')[0].resolution
        ans = 3.46
        self.assertEqual(ans, val)
