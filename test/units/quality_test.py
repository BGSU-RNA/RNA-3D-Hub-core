import unittest
import pytest

from test import StageTest

import pymotifs.units.quality as ntq


class FileHelperTest(StageTest):
    def setUp(self):
        self.helper = ntq.FileHelper()

    def test_can_generate_a_filepath(self):
        val = self.helper('1J5E')
        ans = 'pub/pdb/validation_reports/j5/1j5e/1j5e_validation.xml.gz'
        self.assertEqual(val, ans)


class CoreRsrParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('test/files/validation/4v7w_validation.xml.gz', 'rb') as raw:
            cls.parser = ntq.Parser(raw.read())

    def setUp(self):
        self.parser = self.__class__.parser

    def test_can_generate_a_unit_id(self):
        data = {
            'model': '1',
            'chain': 'A',
            'resname': 'C',
            'resnum': '10',
            'icode': ' '
        }
        val = self.parser._unit_id('1J5E', data)
        ans = {
            'component_id': 'C',
            'chain': 'A',
            'insertion_code': None,
            'component_number': 10,
            'model': 1,
            'pdb': '1J5E',
            'alt_id': None
        }
        self.assertEqual(ans, val)


class MissingRsRParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('test/files/validation/1j5e_validation.xml.gz', 'rb') as raw:
            cls.parser = ntq.Parser(raw.read())

    def setUp(self):
        self.parser = self.__class__.parser

    def test_can_tell_has_no_rsr(self):
        self.assertFalse(self.parser.has_rsr())


class HasRsRParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('test/files/validation/4v7w_validation.xml.gz', 'rb') as raw:
            cls.parser = ntq.Parser(raw.read())

    def setUp(self):
        self.parser = self.__class__.parser

    def test_can_get_tree_from_gz_content(self):
        self.assertTrue(self.parser.root)

    def test_can_detect_has_rsr(self):
        self.assertTrue(self.parser.has_rsr())

    def test_can_generate_nt_level_data(self):
        val = list(self.parser.nts())[0]
        ans = {
            'id': {
                'component_id': 'U',
                'chain': 'AA',
                'insertion_code': None,
                'component_number': 5,
                'model': 1,
                'pdb': '4V7W',
                'alt_id': None,
            },
            'real_space_r': 0.218,
            'z_score': 0.26,
        }
        self.assertEquals(ans, val)


class QueryingTest(StageTest):
    loader_class = ntq.Loader

    def test_knows_if_data_is_missing(self):
        self.assertFalse(self.loader.has_data('0bob'))

    def test_knows_if_data_exists(self):
        self.assertTrue(self.loader.has_data('157D'))


class MappingTest(StageTest):
    loader_class = ntq.Loader

    def test_can_create_complete_mapping(self):
        mapping = self.loader.mapping('1GID')
        self.assertEquals(len(mapping), 350)

    def test_can_map_to_unit_id(self):
        mapping = self.loader.mapping('1GID')
        self.assertEquals(mapping[('A', 146, None, None)], ['1GID|1|A|A|146'])


class MissingDataTest(StageTest):
    loader_class = ntq.Loader

    @classmethod
    def setUpClass(cls):
        with open('test/files/validation/1fjg_validation.xml.gz', 'rb') as raw:
            cls.parser = ntq.Parser(raw.read())
            cls.nts = list(cls.parser.nts())

    def setUp(self):
        self.parser = self.__class__.parser
        self.nts = list(self.__class__.nts)
        super(MissingDataTest, self).setUp()

    def test_parses_all_rna_data(self):
        val = [nt for nt in self.nts if nt['id'].get('chain') == 'A']
        self.assertEquals(len(val), 1603)

    def test_can_map_all_rna_data(self):
        mapping = self.loader.mapping('1JFG')
        val = []
        for nt in self.nts:
            if nt['id']['chain'] == 'A':
                val.append(self.loader.as_quality(nt, mapping))
        self.assertEquals(len(val), 1603)


class AltIdParsingTest(StageTest):
    loader_class = ntq.Loader

    @classmethod
    def setUpClass(cls):
        with open('test/files/validation/1vy4_validation.xml.gz', 'rb') as raw:
            cls.parser = ntq.Parser(raw.read())
            cls.nts = list(cls.parser.nts())

    def test_can_generate_ids_using_alt_ids(self):
        val = [nt for nt in self.nts if nt['id']['alt_id']]
        self.assertTrue(2, len(val))
