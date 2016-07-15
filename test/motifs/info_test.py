import pickle

from test import StageTest

from pymotifs.motifs.info import Loader
from pymotifs.utils import row2dict


class LoaderTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(LoaderTest, self).setUp()
        with open('test/files/motifs/IL.pickle') as raw:
            self.cached = pickle.load(raw)
            self.motifs = self.loader.motifs(self.cached)

    def test_it_computes_all_motifs(self):
        assert len(self.motifs) == 14

    def test_it_assigns_valid_data(self):
        assert row2dict(self.motifs[0]) == {
            'motif_id': 'IL_85752.1',
            'ml_release_id': '0.1',
            'type': 'IL',
            'handle': '85752',
            'version': 1,
            'comment': 'New id, no parents',
        }
