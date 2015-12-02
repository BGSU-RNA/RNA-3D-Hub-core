from test import StageTest

from pymotifs.ife.info import Loader


class QueryTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_data(self):
        self.assertTrue(self.loader.has_data('4V4Q'))

    def test_knows_if_has_no_data(self):
        self.assertFalse(self.loader.has_data('0V4Q'))


class DataTest(StageTest):
    loader_class = Loader

    def test_creates_all_data_for_given_pdb(self):
        val = self.loader.data('4V42')
        # 6 (number of ifes) + 8 (number of chains)
        self.assertEquals(14, len(val))
