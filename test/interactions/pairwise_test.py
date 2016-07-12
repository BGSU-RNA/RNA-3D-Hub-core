import pytest

from test import StageTest
from test import skip_without_matlab

from pymotifs.interactions.pairwise import Loader


class DetectingInteractionTypeTest(StageTest):
    loader_class = Loader

    def test_can_determine_if_bp(self):
        val = self.loader.interaction_type('cHS')
        ans = 'f_lwbp'
        self.assertEquals(val, ans)

    def test_can_determine_if_bph(self):
        val = self.loader.interaction_type('7BPh')
        ans = 'f_bphs'
        self.assertEquals(val, ans)

    def test_can_determine_if_base_ribose(self):
        val = self.loader.interaction_type('0BR')
        ans = 'f_brbs'
        self.assertEquals(val, ans)

    def test_can_determine_if_stack(self):
        val = self.loader.interaction_type('s55')
        ans = 'f_stacks'
        self.assertEquals(val, ans)

    def test_return_none_otherwise(self):
        val = self.loader.interaction_type('cH')
        self.assertTrue(val is None)

    def test_allows_water_as_bp(self):
        val = self.loader.interaction_type('wat')
        ans = 'f_lwbp'
        self.assertEquals(val, ans)

    def test_it_does_not_allow_perp(self):
        self.assertTrue(self.loader.interaction_type('perp') is None)
        self.assertTrue(self.loader.interaction_type('nperp') is None)

    def test_it_does_not_allow_rib(self):
        self.assertTrue(self.loader.interaction_type('rIB') is None)
        self.assertTrue(self.loader.interaction_type('nrIB') is None)


class ParsingACsvTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ParsingACsvTest, self).setUp()
        self.data = self.loader.\
            parse('test/files/interactions/1GID.csv', '1GID')

    def test_can_parse_all_enteris(self):
        val = len(self.data)
        ans = 1513
        self.assertEqual(ans, val)

    def test_can_create_correct_data_structure(self):
        val = self.data[0]
        ans = {
            'unit_id_1': "1GID|1|A|A|104",
            'unit_id_2': "1GID|1|A|A|104",
            'pdb_id': '1GID',
            'f_bphs': '0BPh',
            'f_crossing': 0
        }
        assert val == ans

    @pytest.mark.skip(reason="Not sure of data to use")
    def test_it_merges_entries(self):
        val = self.data[10]  # Not sure what index to use
        ans = {}
        assert val == ans

    def test_it_converts_crossing_to_int(self):
        self.assertEqual(0, self.data[0]['f_crossing'])

    def test_skips_invalid_entries(self):
        ints = [d.get('f_lwbp') for d in self.data]
        self.assertFalse('perp' in ints)


class RunningMatlabTest(StageTest):
    loader_class = Loader

    @skip_without_matlab
    def test_annotates_a_structure(self):
        assert len(self.loader.data('1GID')) == 1513
