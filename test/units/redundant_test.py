from test import StageTest

from pymotifs.units.redundant import RedundantNucleotidesLoader as Loader


class ParsingTest(StageTest):
    loader_class = Loader

    def test_can_parse_a_file(self):
        with open('test/files/redundant-nts.csv', 'rb') as raw:
            data = self.loader._parse(raw, '1GID')
        self.assertEquals(314, len(data))

    def test_can_generate_correct_data(self):
        with open('test/files/redundant-nts.csv', 'rb') as raw:
            data = self.loader._parse(raw, '1GID')
        ans = {
            'unit_id1': "1GID|1|B|G|103",
            'unit_id2': "1GID|1|A|G|103",
            'pdb_id': '1GID'
        }
        self.assertEquals(ans, data[0])


class QueryingTest(StageTest):
    loader_class = Loader

    def test_detects_known_structures(self):
        self.assertTrue(self.loader.has_data('1GID'))

    def test_knows_it_has_no_data(self):
        self.assertFalse(self.loader.has_data('0GID'))
