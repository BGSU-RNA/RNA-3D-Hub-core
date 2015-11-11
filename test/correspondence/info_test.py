from test import StageTest

from pymotifs.correspondence.info import Loader


class LoadingSequencesTest(StageTest):
    loader_class = Loader

    def test_can_load_all_sequenecs(self):
        val = self.loader.exp_seqs(['1G59'])
        ans = [{'id': 1, 'length': 75, 'species': None}]
        self.assertEquals(ans, val)

    def test_can_load_with_synthetic(self):
        val = self.loader.exp_seqs(['1GID'])
        ans = [{'id': 19, 'length': 158, 'species': 32630}]
        self.assertEquals(ans, val)

    def test_can_lookup_several_structures(self):
        val = self.loader.exp_seqs(['1GID', '1G59'])
        ans = [{'id': 1, 'length': 75, 'species': None},
               {'id': 19, 'length': 158, 'species': 32630}]
        self.assertEquals(ans, val)


class HasDataTest(StageTest):
    loader_class = Loader

    def test_can_tell_if_a_pair_was_computed(self):
        val = self.loader.has_data(({'id': 4}, {'id': 5}))
        self.assertTrue(val)

    def test_can_tesll_if_a_pair_was_not_computed(self):
        val = self.loader.has_data(({'id': 1}, {'id': 7}))
        self.assertFalse(val)


class MergingSequencesTest(StageTest):
    loader_class = Loader

    def test_will_create_sets_of_species(self):
        val = self.loader.merge([{'id': 1, 'length': 10, 'species': 10},
                                 {'id': 1, 'length': 10, 'species': None},
                                 {'id': 1, 'length': 10, 'species': 10}])
        ans = [{'id': 1, 'length': 10, 'species': set([10, None])}]
        self.assertEquals(ans, val)

    def test_can_merge_unsorted_chains(self):
        val = self.loader.merge([{'id': 1, 'length': 10, 'species': 10},
                                 {'id': 2, 'length': 1, 'species': 32360},
                                 {'id': 1, 'length': 10, 'species': None},
                                 {'id': 1, 'length': 10, 'species': 10}])
        ans = [{'id': 1, 'length': 10, 'species': set([10, None])},
               {'id': 2, 'length': 1, 'species': set([32360])}]
        self.assertEquals(ans, val)


class ShortSequencesTest(StageTest):
    loader_class = Loader

    def test_will_require_both_equal_if_second_small(self):
        pair = [{'length': 100}, {'length': 30}]
        self.assertFalse(self.loader.valid_length(pair))

    def test_will_pass_if_both_equal(self):
        pair = [{'length': 30}, {'length': 30}]
        self.assertTrue(self.loader.valid_length(pair))

    def test_will_require_both_equal_if_first_small(self):
        pair = [{'length': 10}, {'length': 300}]
        self.assertFalse(self.loader.valid_length(pair))


class LongSequencesTest(StageTest):
    loader_class = Loader

    def test_requires_both_sequences_to_be_above_2000(self):
        pair = [{'length': 1999}, {'length': 2001}]
        self.assertFalse(self.loader.valid_length(pair))

    def test_allows_match_if_both_above_2000(self):
        pair = [{'length': 2001}, {'length': 2002}]
        self.assertTrue(self.loader.valid_length(pair))

    def test_requires_both_sequences_above_36_if_first_below(self):
        pair = [{'length': 35}, {'length': 40}]
        self.assertFalse(self.loader.valid_length(pair))

    def test_requires_both_sequences_above_36_if_second_below(self):
        pair = [{'length': 40}, {'length': 35}]
        self.assertFalse(self.loader.valid_length(pair))

    def test_requires_smallest_atleast_half_larger(self):
        pair = [{'length': 50}, {'length': 102}]
        self.assertFalse(self.loader.valid_length(pair))
        self.assertFalse(self.loader.valid_length(pair[::-1]))

    def test_matches_smallest_atleast_half_larger(self):
        pair = [{'length': 53}, {'length': 102}]
        self.assertTrue(self.loader.valid_length(pair))
        self.assertTrue(self.loader.valid_length(pair[::-1]))


class SpeciesTest(StageTest):
    loader_class = Loader

    def test_requires_species_not_disagree(self):
        pair = [{'species': set([1])}, {'species': set([2])}]
        self.assertFalse(self.loader.valid_species(pair))

    def test_matches_if_one_species_has_None(self):
        pair = [{'species': set([1])}, {'species': set([2, None])}]
        self.assertTrue(self.loader.valid_species(pair))

    def test_matches_if_one_species_has_synenthic(self):
        pair = [{'species': set([1, 32360])}, {'species': set([2])}]
        self.assertTrue(self.loader.valid_species(pair))

    def test_matches_if_same_species(self):
        pair = [{'species': set([1])}, {'species': set([1])}]
        self.assertTrue(self.loader.valid_species(pair))


class ComputingPairsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ComputingPairsTest, self).setUp()
        self.seqs = [
            {'id': 0, 'length': 10, 'species': set([None])},
            {'id': 1, 'length': 10, 'species': set([32360])},
            {'id': 2, 'length': 20, 'species': set([32360])},
            {'id': 3, 'length': 60, 'species': set([30, None])},
            {'id': 4, 'length': 100, 'species': set([562])},
            {'id': 5, 'length': 1000, 'species': set([562])},
            {'id': 6, 'length': 1500, 'species': set([32360, None, 12])},
            {'id': 7, 'length': 1500, 'species': set([12])},
        ]

    def test_can_generate_correct_pairs(self):
        val = self.loader.pairs(self.seqs)
        ans = [
            (self.seqs[0], self.seqs[0]),
            (self.seqs[0], self.seqs[1]),
            (self.seqs[1], self.seqs[1]),
            (self.seqs[2], self.seqs[2]),
            (self.seqs[3], self.seqs[3]),
            (self.seqs[3], self.seqs[4]),
            (self.seqs[4], self.seqs[4]),
            (self.seqs[5], self.seqs[5]),
            (self.seqs[5], self.seqs[6]),
            (self.seqs[6], self.seqs[6]),
            (self.seqs[6], self.seqs[7]),
            (self.seqs[7], self.seqs[7]),
        ]
        self.assertEquals(ans, val)
