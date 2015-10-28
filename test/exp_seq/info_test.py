from test import StageTest

from pymotifs.exp_seq.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_can_get_all_unique_sequences(self):
        self.assertEquals(2, len(self.loader.to_process(['1S72'])))

    def test_knows_if_it_has_data(self):
        self.assertTrue(self.loader.has_data('UUUUCU'))

    def test_knows_if_it_does_not_have_data(self):
        self.assertFalse(self.loader.has_data('UUUUGU'))
