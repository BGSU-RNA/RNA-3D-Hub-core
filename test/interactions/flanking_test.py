import pytest

from test import StageTest

from pymotifs.interactions.flanking import Loader

class ParsingACsvTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ParsingACsvTest, self).setUp()
        self.data = self.loader.\
            parse('test/files/interactions/example-flanking.csv', '2QBG')

    def test_can_parse_all_enteris(self):
        val = len(self.data)
        ans = 5
        self.assertEqual(ans, val)

    def test_can_create_correct_data_structure(self):
        val = self.data[0]
        ans = {
            'unit_id_1': "2QBG|1|B|A|2142",
            'unit_id_2': "2QBG|1|B|U|2149",
            'flanking': 1,
            'pdb_id': '2QBG'
        }
        assert val == ans

    @pytest.mark.skip(reason="Not sure of data to use")
    def test_it_merges_entries(self):
        val = self.data[10]  # Not sure what index to use
        ans = {}
        assert val == ans