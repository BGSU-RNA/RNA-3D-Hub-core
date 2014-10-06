import unittest

import utils as u
import nt_quality as ntq


class FileHelperTest(unittest.TestCase):
    def setUp(self):
        self.helper = ntq.FileHelper()

    def test_can_generate_a_filepath(self):
        val = self.helper('1J5E')
        ans = 'pub/pdb/validation_reports/j5/1j5e/1j5e_validation.xml.gz'
        self.assertEqual(val, ans)


class CoreRsrParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ftp = u.FTPFetchHelper('ftp.wwpdb.org')
        filename = 'pub/pdb/validation_reports/og/3ogy/3ogy_validation.xml.gz'
        cls.parser = ntq.Parser(ftp(filename))

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
        self.assertEqual(val, ans)


class MissingRsRParserTest(unittest.TestCase):
    def setUp(self):
        ftp = u.FTPFetchHelper('ftp.wwpdb.org')
        filename = 'pub/pdb/validation_reports/j5/1j5e/1j5e_validation.xml.gz'
        self.parser = ntq.Parser(ftp(filename))

    def test_can_tell_has_no_rsr(self):
        self.assertFalse(self.parser.has_rsr())


class HasRsRParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ftp = u.FTPFetchHelper('ftp.wwpdb.org')
        filename = 'pub/pdb/validation_reports/og/3ogy/3ogy_validation.xml.gz'
        cls.parser = ntq.Parser(ftp(filename))

    def setUp(self):
        self.parser = self.__class__.parser

    def test_can_get_tree_from_gz_content(self):
        self.assertTrue(self.parser.root)

    def test_can_detect_has_rsr(self):
        self.assertTrue(self.parser.has_rsr())

    def test_can_generate_nt_level_data(self):
        val = list(self.parser.nts())[0]
        ans = {
            'unit_id': '3OGY|1|A|U|5',
            'real_space_r': 0.206
        }
        self.assertEquals(ans, val)
