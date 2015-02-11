from test import StageTest
from pymotifs.chains.best import BestChainsAndModelsLoader


class QueryingTest(StageTest):
    loader_class = BestChainsAndModelsLoader

    def test_knows_if_known_pdb(self):
        self.assertTrue(self.loader.has_data('1S72'))

    def test_knows_if_missingown_pdb(self):
        self.assertFalse(self.loader.has_data('0S72'))
