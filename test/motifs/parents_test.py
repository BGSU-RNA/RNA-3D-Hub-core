import pickle

import pytest

from test import StageTest

from pymotifs.motifs.parents import Loader


class ComputingTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ComputingTest, self).setUp()
        with open('test/files/motifs/IL.pickle') as raw:
            self.cached = pickle.load(raw)
            self.parents = self.loader.parents(self.cached)

    def test_it_can_load_all_parent_data(self):
        assert len(self.parents) == None

    def test_it_creates_correct_data(self):
        assert self.parents[0] == {
            'ml_release_id': '0.1',
            'motif_id': None,
            'parent_ml_release_id': '0.1',
            'parent_motif_id': None,
        }
