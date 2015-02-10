from datetime import datetime
import unittest as ut

from test import StageTest

from pymotifs.pdbs.obsolete import Parser
from pymotifs.pdbs.obsolete import Loader


class ParserTest(ut.TestCase):
    def setUp(self):
        self.parser = Parser()

    def test_can_parse_a_whole(self):
        with open('files/obsolete.dat', 'rb') as raw:
            parsed = self.parser(raw.read())
        self.assertEquals(3226, len(parsed))

    def test_can_extract_correct_data(self):
        with open('files/obsolete.dat', 'rb') as raw:
            parsed = self.parser(raw.read())
        ans = {
            'id': '116L',
            'date': datetime.strptime('31-JUL-94', '%d-%b-%y'),
            'replaced_by': ['216L']
        }
        self.assertEquals(ans, parsed[0])


class ObsoleteTest(StageTest):
    loader_class = Loader

    def test_can_load_the_obsolete_ids(self):
        data = self.loader.data()
        self.assertTrue(len(data) > 0)
