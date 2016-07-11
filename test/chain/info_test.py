from test import StageTest

from pymotifs.chains.info import Loader


class GetIdsTest(StageTest):
    loader_class = Loader

    def test_can_get_known_ids(self):
        reports = [
            {'structureId': '1E4P', 'chainId': 'A'},
            {'structureId': '1DUH', 'chainId': 'A'}
        ]
        reports = self.loader.get_ids(reports)
        val = sorted([report['chain_id'] for report in reports])
        assert len(val) == 2
        types = [type(v) for v in val]
        assert types == [long, long]

    def test_will_not_set_missing_ids(self):
        reports = [{'structureId': '0E4P', 'chainId': 'A'}]
        reports = self.loader.get_ids(reports)
        self.assertTrue('id' not in reports[0])


class GettingDataTest(StageTest):
    loader_class = Loader

    def test_can_get_a_known_pdb(self):
        val = self.loader.data('4V4G')
        self.assertEquals(265, len(val))

    def test_can_process_several_pdbs(self):
        val = self.loader.data(['4V4G', '1S72'])
        self.assertEquals(265 + 31, len(val))


class RenameTest(StageTest):
    loader_class = Loader

    def test_propagates_the_id(self):
        val = self.loader.rename({'chain_id': 'bob', 'pdb_id': 'CX'})
        self.assertEquals('bob', val['chain_id'])

    def test_updates_the_names(self):
        val = self.loader.rename({'pdb_id': 'bob', 'chainId': 'CX'})
        self.assertEquals('CX', val['chain_name'])
        self.assertTrue('chainId' not in val)

    def test_sets_missing_fields_to_none(self):
        val = self.loader.rename({'pdb_id': 'bob', 'chainId': 'CX'})
        self.assertEquals(None, val['sequence'])
