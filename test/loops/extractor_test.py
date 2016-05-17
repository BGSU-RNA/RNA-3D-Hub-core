import os
import shutil

from test import CONFIG
from test import StageTest
from test import skip_without_matlab

from pymotifs.loops.extractor import Loader


class DeterminingComputationTest(StageTest):
    loader_class = Loader

    def test_knows_if_has_pdb_and_type(self):
        self.assertTrue(self.loader.has_data('1S72'))

    def test_knows_if_missing_pdb(self):
        self.assertFalse(self.loader.has_data('0S72'))


class MappingTests(StageTest):
    loader_class = Loader

    def test_can_get_a_correct_mapping_for_HL(self):
        mapping = self.loader._mapping('4V4Q', 'HL')
        self.assertEqual(206, len(mapping))

    def test_can_get_a_correct_mapping_for_IL(self):
        mapping = self.loader._mapping('4V4Q', 'IL')
        self.assertEqual(354, len(mapping))

    def test_can_get_a_correct_mapping_for_J3(self):
        mapping = self.loader._mapping('4V4Q', 'J3')
        self.assertEqual(38, len(mapping))

    def test_gets_empty_mapping_for_missing_structure(self):
        val = self.loader._mapping('0000', 'IL')
        self.assertEquals({}, val)

    def test_has_no_old_style_ids(self):
        mapping = self.loader._mapping('1GID', 'HL')
        seperators = [key[4] for key in mapping.keys()]
        self.assertTrue('_' not in seperators)


class GettingLoopIdsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(GettingLoopIdsTest, self).setUp()
        self.mapping = self.loader._mapping('4V4Q', 'IL')

    def test_generates_a_new_id_for_an_unknown_loop(self):
        val = self.loader._get_loop_id('bob', '4V4Q', 'IL', self.mapping)
        self.assertEquals('IL_4V4Q_355', val)

    def test_adds_units_to_mapping(self):
        val = self.loader._get_loop_id('bob', '4V4Q', 'IL', self.mapping)
        self.assertEquals('IL_4V4Q_355', val)
        self.assertEquals(self.mapping['bob'], 'IL_4V4Q_355')

    def test_uses_old_id_for_known_loop(self):
        nts = (
            '4V4Q|1|AA|A|607,'
            '4V4Q|1|AA|A|608,'
            '4V4Q|1|AA|A|609,'
            '4V4Q|1|AA|A|629,'
            '4V4Q|1|AA|A|630,'
            '4V4Q|1|AA|C|611,'
            '4V4Q|1|AA|C|612,'
            '4V4Q|1|AA|C|631,'
            '4V4Q|1|AA|G|606,'
            '4V4Q|1|AA|G|628,'
            '4V4Q|1|AA|G|633,'
            '4V4Q|1|AA|U|605,'
            '4V4Q|1|AA|U|610,'
            '4V4Q|1|AA|U|632'
        )
        val = self.loader._get_loop_id(nts, '4V4Q', 'IL', self.mapping)
        self.assertEquals('IL_4V4Q_033', val)


class NextNumberTest(StageTest):
    loader_class = Loader

    def test_pads_to_three_for_small_numbers(self):
        val = self.loader._next_loop_number_string(10)
        self.assertEquals('011', val)

    def test_pads_to_six_for_large_numbers(self):
        val = self.loader._next_loop_number_string(1239)
        self.assertEquals('001240', val)

    def test_pading_jumps_at_999(self):
        val = self.loader._next_loop_number_string(999)
        self.assertEquals('001000', val)


class ExtractLoopsTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(ExtractLoopsTest, self).setUp()
        self.loader.save_loops = False
        self.data = self.loader.data('1GID')

    @skip_without_matlab
    def test_can_extract_correct_number_of_loops(self):
        ids = set(d.loop_id for d in self.data)
        names = set(d.loop_name for d in self.data)
        units = set(d.unit_ids for d in self.data)
        self.assertEquals(len(self.data), 22)
        self.assertEqual(len(names), 22)
        self.assertEqual(len(ids), 22)
        self.assertEqual(len(units), 22)


class CreatingFilesTest(StageTest):
    loader_class = Loader
    base = os.path.join(CONFIG['locations']['loops_mat_files'], '1GID')

    def setUp(self):
        super(CreatingFilesTest, self).setUp()
        if os.path.exists(self.base):
            shutil.rmtree(self.base)
        self.data = self.loader.data('1GID')

    @skip_without_matlab
    def test_creates_the_required_files(self):
        def loop(loop_id):
            return str(os.path.join(self.base, loop_id + '.mat'))

        assert os.path.isdir(self.base)
        assert os.path.exists(loop('IL_1GID_001'))
        assert os.path.exists(loop('HL_1GID_001'))

        # 22 loops
        assert len(os.listdir(self.base)) == 22
