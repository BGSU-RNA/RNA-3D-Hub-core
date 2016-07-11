from test import StageTest

from pymotifs.units.coordinates import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_it_knows_if_it_has_data(self):
        assert self.loader.has_data('1GID') is True

    def test_it_knows_if_has_no_data(self):
        assert self.loader.has_data('0GID') is False


class CoordianteTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(CoordianteTest, self).setUp()
        self.data = list(self.loader.data('157D'))

    def test_it_creates_only_atom_level_entries(self):
        for residue in self.data:
            for line in residue['coordinates'].split('\n'):
                assert line.startswith('ATOM')

    def test_it_creates_all_atom_level_entries(self):
        val = self.data[0]['coordinates'].split('\n')
        assert len(val) == 22

    def test_it_creates_entries_for_each_residue(self):
        assert len(self.data) == 24
