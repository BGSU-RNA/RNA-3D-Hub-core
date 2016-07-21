import os

import pytest

from fr3d.unit_ids import decode

from test import StageTest
from test import skip_without_matlab

from pymotifs.loops.positions import Loader


class ParsingTest(StageTest):
    loader_class = Loader

    def test_can_parse_a_file(self):
        positions = self.loader.parse('test/files/loop-positions.csv')
        self.assertEquals(157, len(positions))

    def test_creates_a_correct_data_structure(self):
        positions = self.loader.parse('test/files/loop-positions.csv')
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


class HasDataTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_data(self):
        self.assertTrue(self.loader.has_data('1S72'))

    def test_knows_if_missing_data(self):
        self.assertFalse(self.loader.has_data('0S72'))


class RemovingTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(RemovingTest, self).setUp()

    def test_ignores_unknown_pdb(self):
        self.assertTrue(self.loader.remove('02S7'))

    @pytest.mark.skip()
    def test_deletes_all_positions(self):
        pass


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
            'unit_id': "1GID|1|A|G|149",
            'bulge': 0,
            'flanking': 1,
            'border': 1
        }
        self.assertEquals(ans, self.positions[0])

    def test_will_remove_file_after_loading(self):
        path = './MotifAtlas/Precomputed/1GID/LoopPositions.csv'
        self.assertFalse(os.path.exists(path))


class LoopMappingTest(StageTest):
    loader_class = Loader

    def test_it_builds_a_correct_mapping(self):
        assert self.loader.loop_units_mapping('1DUH') == {
            'IL_1DUH_001': set(['1DUH|1|A|A|39', '1DUH|1|A|A|42', '1DUH|1|A|A|67||||4_555', '1DUH|1|A|A|68||||4_555', '1DUH|1|A|C|40', '1DUH|1|A|C|41', '1DUH|1|A|C|66||||4_555', '1DUH|1|A|G|43', '1DUH|1|A|U|38']),
            'IL_1DUH_002': set(['1DUH|1|A|A|47', '1DUH|1|A|A|59||||4_555', '1DUH|1|A|A|60||||4_555', '1DUH|1|A|A|63||||4_555', '1DUH|1|A|C|46', '1DUH|1|A|C|62||||4_555', '1DUH|1|A|G|48', '1DUH|1|A|G|49', '1DUH|1|A|G|61||||4_555', '1DUH|1|A|G|64||||4_555', '1DUH|1|A|U|45', '1DUH|1|A|U|50']),
            'IL_1DUH_003': set(['1DUH|1|A|A|55', '1DUH|1|A|A|55||||4_555', '1DUH|1|A|A|56', '1DUH|1|A|A|56||||4_555', '1DUH|1|A|C|52', '1DUH|1|A|C|52||||4_555', '1DUH|1|A|G|53', '1DUH|1|A|G|53||||4_555', '1DUH|1|A|G|54', '1DUH|1|A|G|54||||4_555', '1DUH|1|A|G|57', '1DUH|1|A|G|57||||4_555']),
            'IL_1DUH_004': set(['1DUH|1|A|A|47||||4_555', '1DUH|1|A|A|59', '1DUH|1|A|A|60', '1DUH|1|A|A|63', '1DUH|1|A|C|46||||4_555', '1DUH|1|A|C|62', '1DUH|1|A|G|48||||4_555', '1DUH|1|A|G|49||||4_555', '1DUH|1|A|G|61', '1DUH|1|A|G|64', '1DUH|1|A|U|45||||4_555', '1DUH|1|A|U|50||||4_555']),
            'IL_1DUH_005': set(['1DUH|1|A|A|39||||4_555', '1DUH|1|A|A|42||||4_555', '1DUH|1|A|A|67', '1DUH|1|A|A|68', '1DUH|1|A|C|40||||4_555', '1DUH|1|A|C|41||||4_555', '1DUH|1|A|C|66', '1DUH|1|A|G|43||||4_555', '1DUH|1|A|U|38||||4_555']),
        }


class SymmetryOperatorTests(StageTest):
    loader_class = Loader

    def setUp(self):
        super(SymmetryOperatorTests, self).setUp()
        self.positions = self.loader.data('1DUH')

    @pytest.mark.skip(reason="No good data yet")
    def test_it_uses_a_single_symmetry_operator(self):
        units = [pos.unit_id for pos in self.positions]
        sym_ops = set(decode(uid)['symmetry'] for uid in units)
        assert sym_ops == set(['1_555'])
