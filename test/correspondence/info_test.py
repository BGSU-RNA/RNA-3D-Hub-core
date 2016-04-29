import pytest

from test import StageTest

from pymotifs.correspondence.info import Loader


class LoadingSequencesTest(StageTest):
    loader_class = Loader

    def test_can_load_all_sequenecs(self):
        val = self.loader.lookup_sequences('1G59')
        ans = [{'id': 29, 'length': 75, 'species': None}]
        self.assertEquals(ans, val)

    def test_can_load_with_synthetic(self):
        val = self.loader.lookup_sequences('1GID')
        ans = [{'id': 51, 'length': 158, 'species': 32630}]
        self.assertEquals(ans, val)

    @pytest.mark.skip(reason="No data yet")
    def test_it_will_not_load_non_normalized(self):
        self.fail()


class ShortSequencesTest(StageTest):
    loader_class = Loader

    def test_will_require_both_equal_if_second_small(self):
        pair = [{'length': 100}, {'length': 30}]
        self.assertFalse(self.loader.length_match(pair))

    def test_will_pass_if_both_equal(self):
        pair = [{'length': 30}, {'length': 30}]
        self.assertTrue(self.loader.length_match(pair))

    def test_will_require_both_equal_if_first_small(self):
        pair = [{'length': 10}, {'length': 300}]
        self.assertFalse(self.loader.length_match(pair))


class LongSequencesTest(StageTest):
    loader_class = Loader

    def test_requires_both_sequences_to_be_above_2000(self):
        pair = [{'length': 1999}, {'length': 2001}]
        self.assertFalse(self.loader.length_match(pair))

    def test_allows_match_if_both_above_2000(self):
        pair = [{'length': 2001}, {'length': 2002}]
        self.assertTrue(self.loader.length_match(pair))

    def test_requires_both_sequences_above_36_if_first_below(self):
        pair = [{'length': 35}, {'length': 40}]
        self.assertFalse(self.loader.length_match(pair))

    def test_requires_both_sequences_above_36_if_second_below(self):
        pair = [{'length': 40}, {'length': 35}]
        self.assertFalse(self.loader.length_match(pair))

    def test_requires_smallest_atleast_half_larger(self):
        pair = [{'length': 50}, {'length': 102}]
        self.assertFalse(self.loader.length_match(pair))
        self.assertFalse(self.loader.length_match(pair[::-1]))

    def test_matches_smallest_atleast_half_larger(self):
        pair = [{'length': 53}, {'length': 102}]
        self.assertTrue(self.loader.length_match(pair))
        self.assertTrue(self.loader.length_match(pair[::-1]))


class SpeciesTest(StageTest):
    loader_class = Loader

    def test_requires_species_not_disagree(self):
        pair = [{'species': set([1])}, {'species': set([2])}]
        self.assertFalse(self.loader.species_matches(pair))

    def test_does_not_match_if_multiple_disagree(self):
        pair = [{'species': set([1, 3])}, {'species': set([2])}]
        self.assertFalse(self.loader.species_matches(pair))

    def test_matches_if_one_species_has_None(self):
        pair = [{'species': set([1])}, {'species': set([2, None])}]
        self.assertTrue(self.loader.species_matches(pair))

    def test_matches_if_one_species_has_synenthic(self):
        pair = [{'species': set([1, 32360])}, {'species': set([2])}]
        self.assertTrue(self.loader.species_matches(pair))

    def test_matches_if_same_species(self):
        pair = [{'species': set([1])}, {'species': set([1])}]
        self.assertTrue(self.loader.species_matches(pair))

    def test_matches_if_one_in_common(self):
        pair = [{'species': set([1])}, {'species': set([1, 3])}]
        self.assertTrue(self.loader.species_matches(pair))


class IsKnownTest(StageTest):
    loader_class = Loader

    def test_knows_if_a_pair_is_unknown(self):
        pair = [{'id': 1}, {'id': 52}]
        assert self.loader.is_known(pair) is False

    def test_knows_if_a_pair_is_known(self):
        pair = [{'id': 1}, {'id': 2}]
        assert self.loader.is_known(pair) is True


class MatchTest(StageTest):
    loader_class = Loader

    def test_matches_good_length_and_species(self):
        pair = [{'id': -1, 'species': set([1, 32360]), 'length': 53},
                {'id': -2, 'species': set([None]), 'length': 55}]
        assert self.loader.is_match(pair) is True

    def test_does_not_match_if_not_good_species(self):
        pair = [{'id': -1, 'species': set([1]), 'length': 53},
                {'id': -2, 'species': set([3]), 'length': 55}]
        assert self.loader.is_match(pair) is False


class ComputingPairsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ComputingPairsTest, self).setUp()
        pairs = self.loader.pairs(['4V7R', '1GID', '4V88'])
        self.pairs = [(p[0]['id'], p[1]['id']) for p in pairs]

    def test_it_includes_all_self_pairs(self):
        assert (45, 45) in self.pairs
        assert (51, 51) in self.pairs
        assert (50, 50) in self.pairs
        assert (60, 60) in self.pairs
        assert (61, 61) in self.pairs
        assert (62, 62) in self.pairs
        assert (72, 72) in self.pairs

    def test_it_creates_all_possible_pairs(self):
        assert self.pairs == [
            (45, 45), (45, 50), (45, 51), (45, 60), (45, 61), (45, 62), (45, 72),
                      (50, 50), (50, 51), (50, 60), (50, 61), (50, 62), (50, 72),
                                (51, 51), (51, 60), (51, 61), (51, 62), (51, 72),
                                          (60, 60), (60, 61), (60, 62), (60, 72),
                                                    (61, 61), (61, 62), (61, 72),
                                                              (62, 62), (62, 72),
                                                                        (72, 72)
        ]


class ComputingDataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ComputingDataTest, self).setUp()
        self.loader._known = {}
        pairs = self.loader.data(['4V7R', '1GID', '4V88'])
        self.pairs = [(p['exp_seq_id_1'], p['exp_seq_id_2']) for p in pairs]

    def test_includes_comparisions(self):
        assert (45, 50) in self.pairs
        assert (60, 61) in self.pairs
        assert (60, 62) in self.pairs
        assert (61, 62) in self.pairs

    def test_it_includes_self_pairs(self):
        assert (45, 45) in self.pairs
        assert (51, 51) in self.pairs
        assert (50, 50) in self.pairs
        assert (60, 60) in self.pairs
        assert (61, 61) in self.pairs
        assert (62, 62) in self.pairs
        assert (72, 72) in self.pairs

    def test_can_generates_in_correct_ordering(self):
        assert self.pairs == [
            (45, 45),
            (45, 50),
            (50, 50),
            (51, 51),
            (60, 60),
            (60, 61),
            (60, 62),
            (61, 61),
            (61, 62),
            (62, 62),
            (72, 72)
        ]
