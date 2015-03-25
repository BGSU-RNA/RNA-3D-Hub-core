import os

from test import StageTest
from nose import SkipTest

from pymotifs import core
from pymotifs.motifs.base import CsvLoader
from pymotifs.models import MlMutualDiscrepancy


class Example(CsvLoader):
    name = 'example.txt'
    table = MlMutualDiscrepancy
    headers = ['group_name', 'loop_id', 'original_order', 'similarity_order']


class BasicCsvLoaderTests(StageTest):
    loader_class = Example

    def test_knows_current_motif_release(self):
        self.assertEquals('1.15', self.loader.release())

    def test_can_get_base_path(self):
        self.assertEquals('HL_20140614_2131', self.loader.base('HL'))

    def test_can_compute_correct_filename(self):
        raise SkipTest("Need to setup database properly")
        path = os.path.join(os.path.dirname(__file__), "..", "files",
                            "motifs" "IL_20140619_1135", "example.txt")
        self.assertEquals(path, self.loader.filename('IL'))

    def test_will_complain_if_file_is_missing(self):
        self.loader.name = 'missing'
        self.assertRaises(core.InvalidState, self.loader.filename, 'IL')


class ParsingTests(StageTest):
    loader_class = Example

    def setUp(self):
        super(ParsingTests, self).setUp()
        path = os.path.join("test", "files", "motifs", "IL_20120905_0000",
                            "MotifLoopOrder.csv")
        with open(path, 'rb') as raw:
            self.data = self.loader.parse(raw, '1.0')

    def test_can_parse_all_of_a_file(self):
        self.assertEquals(1615, len(self.data))

    def test_can_create_correct_data(self):
        ans = {
            'group_name': 'Group_001',
            'loop_id': 'IL_1ZEV_001',
            'original_order': '1',
            'similarity_order': '74',
            'release_id': '1.0'
        }
        self.assertEquals(ans, self.data[0])
