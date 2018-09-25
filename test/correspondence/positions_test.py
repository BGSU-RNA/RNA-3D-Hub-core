import pytest

from test import StageTest

from pymotifs import core
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.correspondence.positions import Loader


class QueryTest(StageTest):
    loader_class = Loader

    def test_knows_if_correspondence_is_done(self):
        assert self.loader.has_data(1) is True

    def test_knows_if_correspondence_is_not_done(self):
        assert self.loader.has_data(-1) is False

    def test_will_claim_is_not_done_if_length_is_not_null(self):
        old = None
        with self.loader.session() as session:
            corr = session.query(Info).get(1)
            old = corr.length  # Shoud be 4
            assert corr.length is not None
            corr.length = None
            session.merge(corr)

        has = self.loader.has_data(1)

        # Correct data
        with self.loader.session() as session:
            corr = session.query(Info).get(1)
            corr.length = old
            session.merge(corr)

        assert has is False


class LoadingInfoTest(StageTest):
    loader_class = Loader

    def test_it_can_load_info_of_a_correspondence(self):
        assert self.loader.info(1) == (1, 1)

    def test_it_raises_if_unknown(self):
        with pytest.raises(core.InvalidState):
            self.loader.info(-1)


class LoadingSequencesTest(StageTest):
    loader_class = Loader

    def test_it_can_load_an_experimental_sequence(self):
        ans = {'ids': [7388, 7389], 'sequence': 'CA'}
        assert self.loader.sequence(1) == ans

    def test_it_fails_given_invalid_id(self):
        with pytest.raises(core.InvalidState):
            self.loader.sequence(-1)


class CorrelatingTest(StageTest):
    loader_class = Loader

    def test_it_can_correlate_two_sequences(self):
        ref = {'ids': [1, 2, 3], 'sequence': 'AAA'}
        target = {'ids': [10, 11, 12], 'sequence': 'ACA'}
        val = self.loader.correlate(1, ref, target)
        assert val == [
            {
                'exp_seq_position_id_1': 1,
                'exp_seq_position_id_2': 10,
                'correspondence_id': 1,
                'index': 0
            },
            {
                'exp_seq_position_id_1': 2,
                'exp_seq_position_id_2': 11,
                'correspondence_id': 1,
                'index': 1
            },
            {
                'exp_seq_position_id_1': 3,
                'exp_seq_position_id_2': 12,
                'correspondence_id': 1,
                'index': 2
            },
        ]

class ComputingDataTest(StageTest):
    loader_class = Loader

    def test_can_compute_data_in_both_directions(self):
        val = list(self.loader.data(1))
        assert len(val) == 2

    @pytest.mark.skip(reason='No data yet')
    def test_does_not_duplicate_same_alignments(self):
        pass
