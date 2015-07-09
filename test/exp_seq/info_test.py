from test import StageTest

from pymotifs.exp_seq.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_can_get_all_unique_sequences(self):
        self.assertEquals(2, len(self.loader.sequences('1S72')))

    def test_can_generate_all_entries(self):
        self.assertEquals(2, len(self.loader.data('1S72')))
