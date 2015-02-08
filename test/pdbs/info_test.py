from test import StageTest

from pymotifs.pdbs.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_knows_if_data_is_present(self):
        self.assertTrue(self.loader.has_data('2AW7'))

    def test_knows_if_data_is_missing(self):
        self.assertFalse(self.loader.has_data('bob'))

    def test_it_can_create_data(self):
        val = self.loader.data('2AW7')[0].resolution
        ans = 3.46
        self.assertEqual(ans, val)


class RenameTest(StageTest):
    loader_class = Loader

    def test_an_work_if_missing_key(self):
        data = self.loader.rename({'structureId': 'A'})
        self.assertEquals('A', data['id'])

    def test_fills_in_none_for_missing_keys(self):
        data = self.loader.rename({'structureId': 'A'})
        self.assertEquals(None, data['ndb_id'])

    def test_converts_resolution_to_float(self):
        data = self.loader.rename({'resolution': '1'})
        self.assertEquals(1.0, data['resolution'])

    def test_sets_resolution_to_none_if_missing(self):
        data = self.loader.rename({'structureId': '1'})
        self.assertEquals(None, data['resolution'])

    def test_sets_resolution_to_none_if_invalid(self):
        data = self.loader.rename({'resolution': 'A'})
        self.assertEquals(None, data['resolution'])
