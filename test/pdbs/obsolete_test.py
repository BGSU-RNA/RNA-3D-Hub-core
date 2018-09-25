from datetime import datetime
import unittest as ut

from pymotifs.pdbs.obsolete import Parser


class ParserTest(ut.TestCase):
    def setUp(self):
        self.parser = Parser()
        with open('test/files/obsolete.dat', 'rb') as raw:
            self.data = self.parser(raw.read())

    def test_can_parse_a_whole(self):
        self.assertEquals(3226, len(self.data))

    def test_can_extract_correct_data(self):
        ans = {
            'pdb_obsolete_id': '116L',
            'date': datetime.strptime('31-JUL-94', '%d-%b-%y'),
            'replaced_by': '216L'
        }
        self.assertEquals(ans, self.data[0])

    def test_can_replace_empty_with_none(self):
        ans = {
            'pdb_obsolete_id': '179L',
            'date': datetime.strptime("08-JUL-08", '%d-%b-%y'),
            'replaced_by': None
        }
        self.assertEquals(ans, self.data[5])

    def test_can_join_multiple_replacements_to_string(self):
        ans = {
            'pdb_obsolete_id': '1HHB',
            'date': datetime.strptime('18-JUL-84', '%d-%b-%y'),
            'replaced_by': '2HHB,3HHB,4HHB'
        }
        self.assertEquals(ans, self.data[328])
