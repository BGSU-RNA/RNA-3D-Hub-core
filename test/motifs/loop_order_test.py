import pickle

from test import StageTest

from pymotifs.motifs.loop_order import Loader


class ComputingDataTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ComputingDataTest, self).setUp()
        with open('test/files/motifs/v0.1/IL.pickle') as raw:
            self.cached = pickle.load(raw)
            self.ordering = self.loader.ordering(self.cached)

    def test_it_loads_all_data(self):
        assert len(self.ordering) == 81

    def test_it_creates_correct_data(self):
        assert self.ordering[0] == {
            'motif_id': 'IL_85752.1',
            'loop_id': 'IL_1ET4_003',
            'ml_release_id': '0.1',
            'original_order': 1,
            'similarity_order': 1,
        }

    def test_it_creates_valid_data(self):
        assert self.loader.table(**self.ordering[0])
