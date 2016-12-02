import unittest

import pymotifs.quality.utils as ut

import pytest


class ParserTest(unittest.TestCase):
    filename = None

    @classmethod
    def setUpClass(cls):
        with open(cls.filename, 'rb') as raw:
            cls.parser = ut.Parser(raw.read())

    def setUp(self):
        self.parser = self.__class__.parser


class CoreRsrParserTest(ParserTest):
    filename = 'test/files/validation/4v7w_validation.xml.gz'

    def test_can_generate_a_unit_id(self):
        data = {
            'model': '1',
            'chain': 'A',
            'resname': 'C',
            'resnum': '10',
            'icode': ' '
        }
        assert self.parser._unit_id('1J5E', data) == {
            'component_id': 'C',
            'chain': 'A',
            'ins_code': None,
            'number': 10,
            'model': 1,
            'pdb': '1J5E',
            'alt_id': None
        }

    def test_can_get_tree_from_gz_content(self):
        self.assertTrue(self.parser.root)

    def test_can_detect_has_rsr(self):
        self.assertTrue(self.parser.has_rsr())

    def test_can_detect_has_dcc(self):
        assert self.parser.has_dcc() is True

    def test_can_generate_nt_level_data(self):
        assert list(self.parser.nts())[0] == {
            'id': {
                'component_id': 'U',
                'chain': 'AA',
                'ins_code': None,
                'number': 5,
                'model': 1,
                'pdb': '4V7W',
                'alt_id': None,
            },
            'real_space_r': 0.218,
            'real_space_r_z_score': 0.26,
        }

    def test_can_generate_structure_level_data(self):
        assert self.parser.entity() == {
            'pdb_id': '4V7W',
            'percent_rsrz_outliers': 6.76,
            'absolute_percentile_percent_rsrz_outliers': 17.0,
            'relative_percentile_percent_rsrz_outliers': 4.3,
            'clashscore': 40.92,
            'relative_percentile_clashscore': 16.5,
            'absolute_percentile_clashscore': 3.1,
            'percent_rota_outliers': 17.72,
            'absolute_percentile_percent_rota_outliers': 2.5,
            'relative_percentile_percent_rota_outliers': 13.8,
            'md5': 'ad9cd539ce3e7f8c83d1fa706bf3c79a',
        }


class MissingRsRParserTest(ParserTest):
    filename = 'test/files/validation/1j5e_validation.xml.gz'

    def test_can_tell_has_no_rsr(self):
        assert self.parser.has_rsr() is False

    def test_can_tell_has_no_dcc(self):
        assert self.parser.has_dcc() is False


class MissingDataTest(ParserTest):
    filename = 'test/files/validation/1fjg_validation.xml.gz'

    def setUp(self):
        super(MissingDataTest, self).setUp()
        self.parser = self.__class__.parser
        self.nts = list(self.parser.nts())

    def test_parses_all_rna_data(self):
        val = [nt for nt in self.nts if nt['id'].get('chain') == 'A']
        self.assertEquals(len(val), 1603)


class AltIdParsingTest(ParserTest):
    filename = 'test/files/validation/1vy4_validation.xml.gz'

    def setUp(self):
        super(AltIdParsingTest, self).setUp()
        self.parser = self.__class__.parser
        self.nts = list(self.parser.nts())

    def test_can_generate_ids_using_alt_ids(self):
        val = [nt for nt in self.nts if nt['id']['alt_id']]
        self.assertTrue(2, len(val))


class UnusualUnitsTest(ParserTest):
    filename = 'test/files/validation/2uua_validation.xml.gz'

    def setUp(self):
        super(UnusualUnitsTest, self).setUp()
        self.parser = self.__class__.parser
        self.nts = list(self.parser.nts())

    def test_can_generate_unit_ids_for_all_units(self):
        assert len(self.nts) == 4129

    @pytest.mark.skip()
    def test_can_generate_unit_id_for_hard_unit(self):
        data = {
            'real_space_r': 0.32,
            'id': {
                'component_id': 'PAR',
                'chain': 'Z',
                'insertion_code': None,
                'component_number': 1,
                'alt_id': None,
                'model': 1,
                'pdb': '2UUA'
            }
        }
