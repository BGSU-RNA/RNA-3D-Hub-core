from test import StageTest

from pymotifs.correspondence.info import Loader


class LoadingSequencesTest(StageTest):
    loader_class = Loader

    def test_can_load_all_sequenecs(self):
        val = self.loader.lookup_sequences('1G59')
        ans = [
            {'exp_seq_id': 1, 'length': 75, 'species_id': None},
            {'exp_seq_id': 1, 'length': 75, 'species_id': None}
        ]
        self.assertEquals(ans, val)

    def test_can_load_with_synthetic(self):
        val = self.loader.lookup_sequences('1GID')
        ans = [
            {'exp_seq_id': 19, 'length': 158, 'species_id': 32630},
            {'exp_seq_id': 19, 'length': 158, 'species_id': 32630}
        ]
        self.assertEquals(ans, val)


class CheckingPairsTest(StageTest):
    loader_class = Loader

    def test_can_tell_if_a_pair_was_computed(self):
        val = self.loader.is_known({'exp_seq_id1': 4, 'exp_seq_id2': 5})
        self.assertTrue(val)

    def test_can_tesll_if_a_pair_was_not_computed(self):
        val = self.loader.is_known({'exp_seq_id1': 1, 'exp_seq_id2': 7})
        self.assertFalse(val)


class ComputingShortPairsTest(StageTest):
    loader_class = Loader

    def test_will_get_all_same_length_with_no_species(self):
        val = self.loader.pairs({'exp_seq_id': 4, 'length': 8,
                                 'species_id': None})
        ans = [
            {'exp_seq_id1': 4, 'exp_seq_id2': 4},
            {'exp_seq_id1': 4, 'exp_seq_id2': 5}
        ]
        self.assertEquals(ans, val)

    def test_will_get_all_same_length_with_synethic(self):
        val = self.loader.pairs({'exp_seq_id': 4, 'length': 8,
                                 'species_id': 32360})
        ans = [
            {'exp_seq_id1': 4, 'exp_seq_id2': 4},
            {'exp_seq_id1': 4, 'exp_seq_id2': 5}
        ]
        self.assertEquals(ans, val)


class ComputingLongPairsTest(StageTest):
    loader_class = Loader

    def test_can_get_all_similar_length_given_no_species(self):
        val = self.loader.pairs({'exp_seq_id': 33, 'length': 47,
                                 'species_id': None})
        ans = [
            {'exp_seq_id2': 33, 'exp_seq_id1': 1},
            {'exp_seq_id2': 33, 'exp_seq_id1': 6},
            {'exp_seq_id2': 33, 'exp_seq_id1': 9},
            {'exp_seq_id2': 33, 'exp_seq_id1': 11},
            {'exp_seq_id2': 33, 'exp_seq_id1': 14},
            {'exp_seq_id2': 33, 'exp_seq_id1': 24},
            {'exp_seq_id2': 33, 'exp_seq_id1': 26},
            {'exp_seq_id2': 33, 'exp_seq_id1': 30},
            {'exp_seq_id2': 33, 'exp_seq_id1': 32},
            {'exp_seq_id1': 33, 'exp_seq_id2': 33},
            {'exp_seq_id1': 33, 'exp_seq_id2': 46},
            {'exp_seq_id1': 33, 'exp_seq_id2': 48},
            {'exp_seq_id1': 33, 'exp_seq_id2': 49},
        ]
        self.assertEquals(ans, val)

    def test_can_get_all_similar_length_given_synethic(self):
        val = self.loader.pairs({'exp_seq_id': 33, 'length': 47,
                                 'species_id': 32360})
        ans = [
            {'exp_seq_id2': 33, 'exp_seq_id1': 1},
            {'exp_seq_id2': 33, 'exp_seq_id1': 6},
            {'exp_seq_id2': 33, 'exp_seq_id1': 9},
            {'exp_seq_id2': 33, 'exp_seq_id1': 11},
            {'exp_seq_id2': 33, 'exp_seq_id1': 14},
            {'exp_seq_id2': 33, 'exp_seq_id1': 24},
            {'exp_seq_id2': 33, 'exp_seq_id1': 26},
            {'exp_seq_id2': 33, 'exp_seq_id1': 30},
            {'exp_seq_id2': 33, 'exp_seq_id1': 32},
            {'exp_seq_id1': 33, 'exp_seq_id2': 33},
            {'exp_seq_id1': 33, 'exp_seq_id2': 46},
            {'exp_seq_id1': 33, 'exp_seq_id2': 48},
            {'exp_seq_id1': 33, 'exp_seq_id2': 49},
        ]
        self.assertEquals(ans, val)

    def test_can_get_all_similar_length_same_species(self):
        val = self.loader.pairs({'exp_seq_id': 33, 'length': 47,
                                 'species_id': 4932})
        ans = [
            {'exp_seq_id2': 33, 'exp_seq_id1': 1},
            {'exp_seq_id2': 33, 'exp_seq_id1': 6},
            {'exp_seq_id2': 33, 'exp_seq_id1': 9},
            {'exp_seq_id2': 33, 'exp_seq_id1': 11},
            {'exp_seq_id2': 33, 'exp_seq_id1': 14},
            {'exp_seq_id2': 33, 'exp_seq_id1': 24},
            {'exp_seq_id2': 33, 'exp_seq_id1': 32},
            {'exp_seq_id1': 33, 'exp_seq_id2': 33},
        ]
        self.assertEquals(ans, val)
