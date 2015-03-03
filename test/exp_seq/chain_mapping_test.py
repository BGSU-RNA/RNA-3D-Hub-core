from test import StageTest

from pymotifs.exp_seq.chain_mapping import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_gets_correct_chains(self):
        val = self.loader.data('1S72')
        self.assertEquals(2, len(val))
