from test import StageTest

from pymotifs.exp_seq.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_can_get_all_unique_sequences(self):
        self.assertEquals(2128, len(self.loader.possible()))

    def test_can_load_all_known_sequences(self):
        self.assertEquals(0, len(self.loader.known()))

    def test_can_generate_all_entries(self):
        self.assertEquals(2128, len(self.loader.data()))

    def test_can_generate_correct_data(self):
        val = sorted(self.loader.data(), key=lambda s: s.sequence)[0]
        self.assertEquals('AA', val.sequence)
        self.assertEquals(2, val.length)
