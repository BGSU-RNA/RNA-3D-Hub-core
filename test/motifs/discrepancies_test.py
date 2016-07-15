import pickle

import pytest

from test import StageTest

from pymotifs.motifs.discrepancies import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_it_knows_if_version_is_done(self):
        pass

    @pytest.mark.skip()
    def test_it_knows_if_version_is_not_done(self):
        pass


class DiscrepanciesTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(DiscrepanciesTest, self).setUp()
        with open('test/files/motifs/IL.pickle') as raw:
            self.cached = pickle.load(raw)
            self.disc = self.loader.discrepancies(self.cached)

    def test_it_can_load_all_discepancies(self):
        assert len(self.disc) == 2994

    def test_it_creates_correct_data(self):
        assert self.disc[0] == {
            'loop_id_1': 'IL_1F5H_046',
            'loop_id_2': 'IL_1X8W_014',
            'ml_release_id': '0.1',
            'discrepancy': 0.7908,
        }
