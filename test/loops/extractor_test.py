from test import StageTest
from test import skip_without_matlab

from pymotifs.loops.extractor import Loader


class DeterminingComputationTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_pdb_and_type(self):
        self.assertTrue(self.loader.has_data(('1S72', 'IL')))

    def test_knows_if_missing_pdb(self):
        self.assertFalse(self.loader.has_data(('0S72', 'IL')))

    def test_knows_if_missing_type(self):
        self.assertFalse(self.loader.has_data(('1S72', '0L')))

    def test_recomputes_if_missing_loops(self):
        self.assertTrue(self.loader.should_process(('1S72', '0L')))


class GettingLoopIdsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(GettingLoopIdsTest, self).setUp()
        self.mapping = {
            '2AW7_AU_1_A_1076_U_,2AW7_AU_1_A_1077_G_,2AW7_AU_1_A_1078_U_,2AW7_AU_1_A_1079_G_,2AW7_AU_1_A_1080_A_,2AW7_AU_1_A_1081_A_': 'HL_2AW7_024'
        }

    def test_pads_to_three_for_small_numbers(self):
        val = self.loader._next_loop_number_string(10)
        self.assertEquals('011', val)

    def test_pads_to_six_for_large_numbers(self):
        val = self.loader._next_loop_number_string(1239)
        self.assertEquals('001240', val)

    def test_pading_jumps_at_999(self):
        val = self.loader._next_loop_number_string(999)
        self.assertEquals('001000', val)

    def test_generates_a_new_id_for_an_unknown_loop(self):
        val = self.loader._get_loop_id('bob', '2AW7', 'IL', self.mapping, 68)
        self.assertEquals('IL_2AW7_069', val)

    def test_uses_old_id_for_known_loop(self):
        nts = '2AW7_AU_1_A_1076_U_,2AW7_AU_1_A_1077_G_,2AW7_AU_1_A_1078_U_,2AW7_AU_1_A_1079_G_,2AW7_AU_1_A_1080_A_,2AW7_AU_1_A_1081_A_'
        val = self.loader._get_loop_id(nts, '2AW7', 'HL', self.mapping, 10)
        self.assertEquals('HL_2AW7_024', val)


class ExtractLoopsTest(StageTest):
    loader_class = Loader

    @skip_without_matlab
    def test_can_extract_correct_number_of_loops(self):
        val = self.loader.data(('2AW7', 'IL'))
        self.assertEquals(68, len(val))

    @skip_without_matlab
    def test_can_extract_with_reasonable_ids(self):
        data = self.loader.data(('4V4Q', 'IL'))
        names = [data.name for name in data]
        self.assertTrue(len(set(names)) > 69)
