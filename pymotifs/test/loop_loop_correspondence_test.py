import unittest as ut

from models import Session
import loop_loop_correspondence as ll


class StructureUtilTest(ut.TestCase):
    def setUp(self):
        self.util = ll.StructureUtil(Session)

    def test_gets_all_loops(self):
        val = len(self.util.loops('2AW7'))
        ans = 109
        self.assertEquals(ans, val)

    def test_loads_loop_data(self):
        val = self.util.loops('2AW7')[0]
        ans = {
            'id': 'HL_2AW7_001',
            'nts': ['2AW7_AU_1_A_12_U_', '2AW7_AU_1_A_13_U_',
                    '2AW7_AU_1_A_14_U_', '2AW7_AU_1_A_15_G_',
                    '2AW7_AU_1_A_16_A_', '2AW7_AU_1_A_17_U_',
                    '2AW7_AU_1_A_18_C_', '2AW7_AU_1_A_19_A_',
                    '2AW7_AU_1_A_20_U_', '2AW7_AU_1_A_21_G_',
                    '2AW7_AU_1_A_22_G_']
        }
        self.assertEquals(ans, val)

    def test_gets_reference_structures(self):
        val = self.util.reference('3J5F')
        ans = ['4KJA', '2AW7']
        self.assertEquals(ans, val)


class StructureUtilMappingTest(ut.TestCase):
    def setUp(self):
        self.util = ll.StructureUtil(Session)
        self.mapping = self.util.mapping('4KJA', '3J55')

    def test_gives_empty_mapping_for_invalid_pair(self):
        val = self.util.mapping('1J5E', 'bob')
        self.assertEquals({}, val)

    def test_gets_mapping_from_ref_to_given(self):
        val = self.mapping['4KJA_AU_1_A_13_U_']
        ans = '3J55_AU_1_A_13_U_'
        self.assertEqual(ans, val)

    def test_gets_mapping_from_given_to_ref(self):
        val = self.mapping['3J55_AU_1_A_24_U_']
        ans = '4KJA_AU_1_A_24_U_'
        self.assertEqual(ans, val)
