import unittest

from models import Session

import correspondence.nts as ntnt
import utils as ut

VALID_CORRELATIONS = '''[
{"unit1": "A", "unit2": "1"},
{"unit1": "C", "unit2": "2"}
{"unit1": "-", "unit2": "2"}
{"unit1": "C", "unit2": "-"}
]'''


class CorrelationResponseParserTest(unittest.TestCase):
    def setUp(self):
        self.parser = ntnt.Parser()

    def test_can_parse_valid_response(self):
        val = self.parser(VALID_CORRELATIONS)
        ans = [{'unit1_id': 'A', 'unit2_id': '1'},
               {'unit1_id': 'C', 'unit2_id': '2'}]
        self.assertEqual(ans, val)

    def test_can_add_additional_data(self):
        self.parser.additional = {'bob': 'jones'}
        val = self.parser(VALID_CORRELATIONS)
        ans = [{'unit1_id': 'A', 'unit2_id': '1', 'bob': 'jones'},
               {'unit1_id': 'C', 'unit2_id': '2', 'bob': 'jones'}]
        self.assertEqual(ans, val)

    def test_fails_with_empty_response(self):
        self.assertRaises(ut.EmptyResponse, self.parser, '')


class StructureUtilTest(unittest.TestCase):
    def setUp(self):
        self.util = ntnt.StructureUtil(Session)

    def test_can_find_representative(self):
        val = self.util.representative('1J5E')
        self.assertEquals(set(['2VQE', '1FJG']), val)

    def test_does_not_give_self_as_representataive(self):
        val = self.util.representative('2VQE')
        self.assertEquals(set(['1FJG']), val)

    def test_can_get_longest_chain(self):
        val = self.util.longest_chain('2AW7')
        self.assertEquals('A', val)

    def test_can_get_sequences_to_correlate(self):
        ans = ['UGAAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCUGGGAAACUGCCUGAUGGAGGGGGAUAACUACUGGAAACGGUAGCUAAUACCGCAUAACGUCGCAAGACCAAAGAGGGGGACCUUCGGGCCUCUUGCCAUCGGAUGUGCCCAGAUGGGAUUAGCUAGUAGGUGGGGUAACGGCUCACCUAGGCGACGAUCCCUAGCUGGUCUGAGAGGAUGACCAGCCACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGCGCAAGCCUGAUGCAGCCAUGCCGCGUGUAUGAAGAAGGCCUUCGGGUUGUAAAGUACUUUCAGCGGGGAGGAAGGGAGUAAAGUUAAUACCUUUGCUCAUUGACGUUACCCGCAGAAGAAGCACCGGCUAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGUGCAAGCGUUAAUCGGAAUUACUGGGCGUAAAGCGCACGCAGGCGGUUUGUUAAGUCAGAUGUGAAAUCCCCGGGCUCAACCUGGGAACUGCAUCUGAUACUGGCAAGCUUGAGUCUCGUAGAGGGGGGUAGAAUUCCAGGUGUAGCGGUGAAAUGCGUAGAGAUCUGGAGGAAUACCGGUGGCGAAGGCGGCCCCCUGGACGAAGACUGACGCUCAGGUGCGAAAGCGUGGGGAGCAAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGAUGUCGACUUGGAGGUUGUGCCCUUGAGGCGUGGCUUCCGGAGCUAACGCGUUAAGUCGACCGCCUGGGGAGUACGGCCGCAAGGUUAAAACUCAAAUGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAUGCAACGCGAAGAACCUUACCUGGUCUUGACAUCCACGGAAGUUUUCAGAGAUGAGAAUGUGCCUUCGGGAACCGUGAGACAGGUGCUGCAUGGCUGUCGUCAGCUCGUGUUGUGAAAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUAUCCUUUGUUGCCAGCGGUCCGGCCGGGAACUCAAAGGAGACUGCCAGUGAUAAACUGGAGGAAGGUGGGGAUGACGUCAAGUCAUCAUGGCCCUUACGACCAGGGCUACACACGUGCUACAAUGGCGCAUACAAAGAGAAGCGACCUCGCGAGAGCAAGCGGACCUCAUAAAGUGCGUCGUAGUCCGGAUUGGAGUCUGCAACUCGACUCCAUGAAGUCGGAAUCGCUAGUAAUCGUGGAUCAGAAUGCCACGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGGUAGCUUAACCUUCGGGAGGGCGCUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAGGUAACCGUAGGGGAACCUGCGGUUGGAUCA']
        val = self.util.polymer_sequences('2AW7', 'A')
        self.assertEqual(ans, val)

    def test_it_can_get_all_ids(self):
        val = self.util.unit_ids('2AW7', 'A')
        self.assertEqual(1530, len(val))


class NtNtCorrespondencesTest(unittest.TestCase):
    def setUp(self):
        self.corr = ntnt.Loader({}, Session)

    def test_can_check_if_has_correspondence(self):
        self.assertFalse(self.corr.has_correspondence('1J5E', 'bob'))
