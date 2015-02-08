from datetime import datetime
import unittest as ut

from test import StageTest

from pymotifs.pdbs.obsolete import Parser
from pymotifs.pdbs.obsolete import Loader


class ParserTest(ut.TestCase):
    def setUp(self):
        self.parser = Parser()

    def test_can_parse_a_single_line(self):
        val = self.parser(['OBSLTE    26-SEP-06 2H33     2JM5 2OWI'])
        ans = [{
            'id': '2H33',
            'date': datetime.strptime('26-SEP-06', '%d-%b-%y'),
            'replaced_by': ['2JM5', '2OWI']
        }]
        self.assertEquals(ans, val)

    def test_can_parse_a_whole(self):
        with open('files/obsolete.dat', 'rb') as raw:
            parsed = self.parser(raw)
        self.assertEquals(3226, len(parsed))


class ObsoleteTest(StageTest):
    loader_class = Loader

    def test_can_load_the_obsolete_ids(self):
        data = self.loader.data()
        self.assertTrue(len(data) > 0)
