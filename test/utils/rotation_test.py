from test import StageTest

from pymotifs.units.rotation import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_it_knows_if_has_not_been_done(self):
        assert self.loader.has_data('0FGJ') is False

    def test_it_knows_if_has_been_done(self):
        assert self.loader.has_data('1IBK') is True
