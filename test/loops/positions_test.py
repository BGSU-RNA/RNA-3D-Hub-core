import os

from test import StageTest
from test import skip_without_matlab

from pymotifs.loops.positions import Loader


class ParsingTest(StageTest):
    loader_class = Loader

    def test_can_parse_a_file(self):
        positions = self.loader.parse('files/loop-positions.csv')
        self.assertEquals(157, len(positions))

    def test_creates_a_correct_data_structure(self):
        positions = self.loader.parse('files/loop-positions.csv')
        val = positions[0]
        ans = {
            'loop_id': "HL_1GID_001",
            'position': 1,
            'unit_id': "1GID|1|B|G|234",
            'bulge': 0,
            'flanking': 1,
            'border': 1
        }
        self.assertEquals(ans, val)


@skip_without_matlab
class MatlabTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(MatlabTest, self).setUp()
        self.positions = self.loader.annotations('1GID')

    def test_generates_all_data(self):
        self.assertEquals(194, len(self.positions))

    def test_creates_valid_data_structures(self):
        ans = {
            'loop_id': "HL_1GID_001",
            'position': 1,
            'unit_id': "1GID|1|B|G|234",
            'bulge': 0,
            'flanking': 1,
            'border': 1
        }
        self.assertEquals(ans, self.positions[0])

    def test_will_remove_file_after_loading(self):
        path = './MotifAtlas/Precomputed/1GID/LoopPositions.csv'
        self.assertFalse(os.path.exists(path))
