import pytest

from test import StageTest

from pymotifs import models as mod
from pymotifs.constants import SYNTHENIC_SPECIES_ID
from pymotifs.correspondence.info import Loader


class LoadingSequencesTest(StageTest):
    loader_class = Loader

    def test_can_load_all_sequenecs(self):
        val = self.loader.lookup_sequences('1G59')
        assert val == [{'id': 24, 'length': 75, 'species': None}]

    def test_can_load_with_synthetic(self):
        val = self.loader.lookup_sequences('1GID')
        assert val == [{'id': 40, 'length': 158, 'species': SYNTHENIC_SPECIES_ID}]

    @pytest.mark.skip(reason="No data yet")
    def test_it_will_not_load_non_normalized(self):
        self.fail()


class ShortSequencesTest(StageTest):
    loader_class = Loader

    def test_will_require_both_equal_if_second_small(self):
        pair = [{'length': 100}, {'length': 30}]
        assert self.loader.length_match(pair) is False

    def test_will_pass_if_both_equal(self):
        pair = [{'length': 30}, {'length': 30}]
        assert self.loader.length_match(pair) is True

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

    def test_requires_both_sequences_above_19_if_first_below(self):
        pair = [{'length': 18}, {'length': 20}]
        self.assertFalse(self.loader.length_match(pair))

    def test_requires_both_sequences_above_19_if_second_below(self):
        pair = [{'length': 20}, {'length': 18}]
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
        pair = [{'species': set([1, SYNTHENIC_SPECIES_ID])},
                {'species': set([2])}]
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
        assert self.loader.is_known([{'id': 1}, {'id': 52}]) is False

    def test_knows_if_a_pair_is_known(self):
        assert self.loader.is_known([{'id': 1}, {'id': 68}]) is True


class MatchTest(StageTest):
    loader_class = Loader

    def test_matches_good_length_and_species(self):
        pair = [{'id': -1, 'species': set([1, SYNTHENIC_SPECIES_ID]), 'length': 53},
                {'id': -2, 'species': set([None]), 'length': 55}]
        assert self.loader.is_match(pair) is True

    def test_does_not_match_if_not_good_species(self):
        pair = [{'id': -1, 'species': set([1]), 'length': 53},
                {'id': -2, 'species': set([3]), 'length': 55}]
        assert self.loader.is_match(pair) is False


class ComputingDataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ComputingDataTest, self).setUp()
        self.loader._known = {}
        pairs = self.loader.data(['4V7R', '1GID', '4V88'])
        self.pairs = [(p['exp_seq_id_1'], p['exp_seq_id_2']) for p in pairs]

    def exp_seq_id(self, pdb, chain):
        with self.loader.session() as session:
            return session.query(mod.ExpSeqPdb).\
                filter_by(pdb_id=pdb, chain_name=chain).\
                one().\
                exp_seq_id

    def exp_pair(self, pdb1, chain1, pdb2, chain2):
        return tuple(sorted([self.exp_seq_id(pdb1, chain1),
                             self.exp_seq_id(pdb2, chain2)]))

    def self_pair(self, pdb, chain):
        return (self.exp_seq_id(pdb, chain), self.exp_seq_id(pdb, chain))

    def test_includes_comparisions(self):
        assert self.exp_pair('4V7R', 'A1', '4V88', 'A2') in self.pairs
        assert self.exp_pair('4V7R', 'B1', '4V88', 'A1') in self.pairs
        assert self.exp_pair('4V7R', 'B2', '4V88', 'A3') in self.pairs
        assert self.exp_pair('4V7R', 'B3', '4V88', 'A4') in self.pairs
        assert self.exp_pair('4V7R', 'C1', '4V88', 'A6') in self.pairs
        assert self.exp_pair('1GID', 'A', '4V88', 'A3') in self.pairs
        assert self.exp_pair('1GID', 'A', '4V88', 'A4') in self.pairs
        assert self.exp_pair('1GID', 'A', '4V7R', 'B2') in self.pairs
        assert self.exp_pair('1GID', 'A', '4V7R', 'B3') in self.pairs
        assert self.exp_pair('4V88', 'A2', '4V88', 'A6') in self.pairs
        assert self.exp_pair('4V7R', 'B2', '4V7R', 'B3') in self.pairs

    def test_it_includes_self_pairs(self):
        assert self.self_pair('4V7R', 'A1') in self.pairs
        assert self.self_pair('4V7R', 'B1') in self.pairs
        assert self.self_pair('4V7R', 'B2') in self.pairs
        assert self.self_pair('4V7R', 'B3') in self.pairs
        assert self.self_pair('1GID', 'A') in self.pairs
        assert self.self_pair('4V88', 'A2') in self.pairs
        assert self.self_pair('4V88', 'A1') in self.pairs
        assert self.self_pair('4V88', 'A2') in self.pairs
        assert self.self_pair('4V88', 'A3') in self.pairs
        assert self.self_pair('4V88', 'A4') in self.pairs
        assert self.self_pair('4V88', 'A6') in self.pairs

    def test_can_generates_in_correct_ordering(self):
        ans = sorted(set([
            self.self_pair('4V7R', 'A1'),
            self.self_pair('4V7R', 'C1'),
            self.self_pair('4V7R', 'B1'),
            self.self_pair('4V7R', 'B2'),
            self.self_pair('4V7R', 'B3'),
            self.self_pair('1GID', 'A'),
            self.self_pair('4V88', 'A2'),
            self.self_pair('4V88', 'A1'),
            self.self_pair('4V88', 'A2'),
            self.self_pair('4V88', 'A3'),
            self.self_pair('4V88', 'A4'),
            self.self_pair('4V88', 'A5'),
            self.self_pair('4V88', 'A6'),
            self.exp_pair('4V7R', 'A1', '4V88', 'A2'),
            self.exp_pair('4V7R', 'C1', '4V88', 'A6'),
            self.exp_pair('4V7R', 'B1', '4V88', 'A1'),
            self.exp_pair('4V7R', 'B2', '4V88', 'A3'),
            self.exp_pair('4V7R', 'B3', '4V88', 'A4'),
            self.exp_pair('1GID', 'A', '4V88', 'A3'),
            self.exp_pair('1GID', 'A', '4V7R', 'B2'),
            self.exp_pair('1GID', 'A', '4V7R', 'B3'),
            self.exp_pair('1GID', 'A', '4V88', 'A4'),
            self.exp_pair('4V88', 'A2', '4V88', 'A6'),
            self.exp_pair('4V7R', 'B2', '4V7R', 'B3'),
        ]))
        print(self.pairs)
        print(ans)
        assert self.pairs == ans
