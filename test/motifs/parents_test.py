import pickle

import pytest

from test import StageTest

from pymotifs import core
from pymotifs.motifs.parents import Loader


class BaseTest(StageTest):
    loader_class = Loader

    def data(self, release):
        with open('test/files/motifs/v%s/IL.pickle' % release, 'rb') as raw:
            cached = pickle.load(raw)
            return self.loader.parents(cached)


class NoParentsTest(BaseTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_it_knows_when_has_no_parents(self):
        assert self.loader.no_parents(self.data('0.1')) is True

    @pytest.mark.skip()
    def test_raises_skip_when_no_parents(self):
        with pytest.raises(core.Skip):
            self.loader.data(self.data('0.1'))

    @pytest.mark.skip()
    def test_gives_empty_list_for_no_parents(self):
        assert self.loader.data(self.data('0.1')) == []


class ParentsTest(BaseTest):
    @pytest.mark.skip()
    def test_it_creates_correct_data(self):
        assert self.parents[0] == {
            'ml_release_id': '0.1',
            'motif_id': None,
            'parent_ml_release_id': '0.1',
            'parent_motif_id': None,
        }
