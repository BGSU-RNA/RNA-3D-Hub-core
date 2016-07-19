import pickle

import pytest

from test import StageTest

from pymotifs.motifs.loop_positions import Loader


class ComputingDataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ComputingDataTest, self).setUp()
        with open('test/files/motifs/IL.pickle') as raw:
            self.cached = pickle.load(raw)
            self.positions = self.loader.positions(self.cached)

    def test_it_loads_all_data(self):
        assert len(self.positions) == 527

    def test_it_creates_correct_data(self):
        assert self.positions[0] == {
            'motif_id': 'IL_85752.1',
            'loop_id': 'IL_1ET4_003',
            'ml_release_id': '0.1',
            'position': 1,
            'unit_id': '1ET4|1|B|C|311',
        }

    def test_it_creates_valid_data(self):
        assert self.loader.table(**self.positions[0])
