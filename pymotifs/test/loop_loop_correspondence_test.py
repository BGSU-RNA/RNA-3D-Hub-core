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


class StructureCoverageTest(ut.TestCase):
    def setUp(self):
        self.util = ll.StructureUtil(Session)
        self.loader = ll.Loader(Session)
        self.mapping = self.util.mapping('4KJA', '3J55')

    def test_fails_if_given_empty_first_loop(self):
        loop = {'id': 'A', 'nts': ['4KJA_AU_1_A_18_C_']}
        val = self.loader.coverage({}, loop, self.mapping)
        self.assertTrue(val is None)

    def test_fails_if_given_empty_second_loop(self):
        loop = {'id': 'A', 'nts': ['4KJA_AU_1_A_18_C_']}
        val = self.loader.coverage(loop, {}, self.mapping)
        self.assertTrue(val is None)

    def test_fails_if_first_nts_empty(self):
        loop = {'id': 'A', 'nts': ['4KJA_AU_1_A_18_C_']}
        empty = {'id': None, 'nts': []}
        val = self.loader.coverage(empty, loop, self.mapping)
        self.assertTrue(val is None)

    def test_fails_if_second_nts_empty(self):
        loop = {'id': 'A', 'nts': ['4KJA_AU_1_A_18_C_']}
        empty = {'id': None, 'nts': []}
        val = self.loader.coverage(loop, empty, self.mapping)
        self.assertTrue(val is None)

    def test_fails_if_one_does_not_map(self):
        loop = {
            'id': 'A',
            'nts': ['4KJA_AU_1_Q_18_C_', '4KJA_AU_1_A_18_C_']
        }
        val = self.loader.coverage(loop, loop, self.mapping)
        self.assertTrue(val is None)

    def test_fails_if_all_do_not_map(self):
        loop = {'id': 'A', 'nts': ['4KJA_AU_1_Q_18_C_']}
        val = self.loader.coverage(loop, loop, self.mapping)
        self.assertTrue(val is None)

    def test_gets_exact_coverage(self):
        pass
        # loop1 = {'id': 'A', 'nts': ['4KJA_AU_1_Q_18_C_']}
        # loop2 = {'id': 'A', 'nts': ['']}
        # val = self.loader.coverage(loop1, loop2, self.mapping)
        # self.assertEqual('exact', val)

    def test_gets_partial(self):
        pass

    def test_gets_parital_if_one_in_common(self):
        pass

    def test_finds_contained(self):
        pass

    def test_finds_enclose(self):
        pass


class CompareTests(ut.TestCase):
    def setUp(self):
        self.loader = ll.Loader(Session)

    def test_process_all_loops_in_both_structures(self):
        pass
        # val = len(list(self.loader.compare('2AW7', '4KJA')))
        # ans = None
        # self.assertEqual(ans, val)

    def test_assigns_unique_to_unique_loop(self):
        pass

    def test_generates_valid_objects(self):
        pass
