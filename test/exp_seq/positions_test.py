from test import StageTest

from pymotifs.exp_seq.positions import Loader


class GeneratingDataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(GeneratingDataTest, self).setUp()
        self.data = self.loader.positions(10, 'ACGI', 'ACGG')

    def test_can_find_all_positions(self):
        self.assertEquals(4, len(self.data))

    def test_can_generate_valid_data(self):
        ans = {
            'exp_seq_id': 10,
            'unit': 'I',
            'normalized_unit': 'G',
            'index': 3,
        }
        self.assertEquals(ans, self.data[3])


class GettingSequencesTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(GettingSequencesTest, self).setUp()
        self.data = sorted(self.loader.to_process(['1S72']))

    def test_finds_all_sequences(self):
        self.assertEquals(2, len(self.data))
