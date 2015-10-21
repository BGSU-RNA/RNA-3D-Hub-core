from test import StageTest

from pymotifs.pdbs.info import Loader


class InfoTest(StageTest):
    loader_class = Loader

    def test_it_can_create_data(self):
        val = self.loader.data(['2AW7'])[0].resolution
        ans = 3.46
        self.assertEqual(ans, val)

    def test_it_can_load_several_structures(self):
        val = self.loader.data(['2AW7', '1GID', '1S72'])
        self.assertEquals(3, len(val))

    def test_it_gets_the_data_correctly(self):
        data = self.loader.data(['2AW7', '1GID', '1S72'])
        val = sorted([d.resolution for d in data])
        self.assertEquals([2.4, 2.5, 3.46], val)


class RenameTest(StageTest):
    loader_class = Loader

    def test_an_work_if_missing_key(self):
        data = self.loader.rename({'structureId': 'A'})
        self.assertEquals('A', data['pdb_id'])

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

    def test_sets_resolution_to_none_if_empty(self):
        data = self.loader.rename({'resolution': ''})
        self.assertEquals(None, data['resolution'])
