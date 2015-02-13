import unittest
from nose import SkipTest

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
        with open('files/validation/4v7w_validation.xml.gz', 'rb') as raw:
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
        ans = '1J5E|1|A|C|10'
        self.assertEqual(ans, val)


class MissingRsRParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/validation/1j5e_validation.xml.gz', 'rb') as raw:
            cls.parser = ntq.Parser(raw.read())

    def setUp(self):
        self.parser = self.__class__.parser

    def test_can_tell_has_no_rsr(self):
        self.assertFalse(self.parser.has_rsr())


class HasRsRParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/validation/4v7w_validation.xml.gz', 'rb') as raw:
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
            'id': '4V7W|1|AA|U|5',
            'real_space_r': 0.218
        }
        self.assertEquals(ans, val)


class QueryingTest(StageTest):
    loader_class = ntq.Loader

    def test_knows_if_data_is_missing(self):
        self.assertFalse(self.loader.has_data('0bob'))

    def test_knows_if_data_exists(self):
        raise SkipTest("No good data")
