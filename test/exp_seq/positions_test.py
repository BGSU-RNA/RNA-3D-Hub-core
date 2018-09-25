from test import StageTest

from pymotifs.exp_seq.positions import Loader


class GeneratingDataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(GeneratingDataTest, self).setUp()
        self.data = self.loader.positions(10, 'ACGI')

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


class NotNormalizedSequencetest(StageTest):
    loader_class = Loader

    def test_can_deal_with_nonnormalized_sequence(self):
        val = self.loader.positions(1, 'AAF')
        assert val == [
            {'exp_seq_id': 1, 'unit': 'A', 'normalized_unit': 'A', 'index': 0},
            {'exp_seq_id': 1, 'unit': 'A', 'normalized_unit': 'A', 'index': 1},
            {'exp_seq_id': 1, 'unit': 'F', 'normalized_unit': None, 'index': 2},
        ]


class RealDataTest(StageTest):
    loader_class = Loader

    def test_knows_if_sequence_is_processed(self):
        assert self.loader.has_data(1) is True

    def test_knows_if_sequence_is_not_processed(self):
        assert self.loader.has_data(-1) is False

    def test_can_correctly_get_exp_seq_ids(self):
        assert self.loader.to_process(['1GID']) == [40]

    def test_it_can_load_a_sequence(self):
        assert self.loader.sequence(1) == 'CA'

    def test_can_create_positions(self):
        assert self.loader.data(1) == [
            {'exp_seq_id': 1, 'unit': 'C', 'normalized_unit': 'C', 'index': 0},
            {'exp_seq_id': 1, 'unit': 'A', 'normalized_unit': 'A', 'index': 1},
        ]
