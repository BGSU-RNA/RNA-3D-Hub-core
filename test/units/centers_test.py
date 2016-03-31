import pytest

from test import StageTest
from test import CifStageTest

from pymotifs.units.centers import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_it_has_no_data(self):
        assert self.loader.has_data('024D') is False

    @pytest.mark.xfail()
    def test_knows_it_has_data(self):
        assert self.loader.has_data('124D') is True

    @pytest.mark.skip()
    def test_can_delete_data(self):
        pass


class LoadingCenters(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/124D.cif'

    def test_it_can_compute_centers(self):
        val = list(self.loader.data(self.structure))
        assert len(val) > 0
