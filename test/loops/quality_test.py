import pytest

from pymotifs import core

from test import StageTest

from pymotifs.loops.quality import Loader


class Base(StageTest):
    def loop(self, loop_id):
        pdb = loop_id.split('_')[1]
        for loop in self.loader.loops(pdb):
            if loop['id'] == loop_id:
                return loop
        self.fail("Could not find loop %s" % loop_id)


class QueryingTest(Base):
    loader_class = Loader

    def test_can_detect_if_has_data(self):
        self.assertTrue(self.loader.has_data('1GID'))

    def test_can_detect_if_has_no_data(self):
        self.assertFalse(self.loader.has_data('0GID'))


class ComplementarySequenceTest(Base):
    loader_class = Loader

    def test_can_detect_if_has_non_cWW_pairs(self):
        loop = self.loop('HL_1GID_001')
        assert self.loader.has_no_non_cWW(loop) is False

    def test_can_find_if_has_no_non_cWW_pairs(self):
        loop = self.loop('IL_1GRZ_007')
        assert self.loader.has_no_non_cWW(loop) is True

    def test_can_find_if_complemenatry_sequence(self):
        valid = {'units': ['A', 'C', 'C', 'G']}
        invalid1 = {'units': ['A', 'C', 'G', 'U']}
        invalid2 = {'units': ['A', 'C', 'A', 'G', 'U']}
        assert self.loader.complementary_sequence(valid) is False
        assert self.loader.complementary_sequence(invalid1) is True
        assert self.loader.complementary_sequence(invalid2) is False

    def test_knows_if_not_complementary(self):
        loop = self.loop('IL_3CPW_080')
        assert self.loader.has_no_non_cWW(loop) is True
        assert self.loader.complementary_sequence(loop) is False
        assert self.loader.is_complementary(loop) is False

    def test_handles_multiple_parts(self):
        loop = self.loop('IL_1S72_068')
        assert self.loader.has_no_non_cWW(loop) is True
        assert self.loader.complementary_sequence(loop) is False
        assert self.loader.is_complementary(loop) is False


class ChainNumberValidationTests(Base):
    loader_class = Loader

    def bad_chain_number(self, *args, **kwargs):
        return self.loader.bad_chain_number(*args, **kwargs)

    def test_knows_if_good_IL(self):
        loop1 = {'type': 'IL', 'chains': set('A')}
        loop2 = {'type': 'IL', 'chains': set('AB')}
        assert self.bad_chain_number(loop1) is False
        assert self.bad_chain_number(loop2) is False

    def test_knows_if_has_bad_IL(self):
        loop1 = {'type': 'IL', 'chains': set('ABC')}
        loop2 = {'type': 'IL', 'chains': set('')}
        assert self.bad_chain_number(loop1) is True
        assert self.bad_chain_number(loop2) is True

    def test_knows_if_good_HL(self):
        loop1 = {'type': 'HL', 'chains': set('A')}
        assert self.bad_chain_number(loop1) is False

    def test_knows_if_bad_HL(self):
        loop1 = {'type': 'HL', 'chains': set('ABC')}
        loop2 = {'type': 'HL', 'chains': set('BC')}
        loop3 = {'type': 'HL', 'chains': set('')}
        assert self.bad_chain_number(loop1) is True
        assert self.bad_chain_number(loop2) is True
        assert self.bad_chain_number(loop3) is True

    def test_knows_if_good_J3(self):
        loop1 = {'type': 'J3', 'chains': set('A')}
        loop2 = {'type': 'J3', 'chains': set('AB')}
        loop3 = {'type': 'J3', 'chains': set('ABC')}
        assert self.bad_chain_number(loop1) is False
        assert self.bad_chain_number(loop2) is False
        assert self.bad_chain_number(loop3) is False

    def test_knows_if_bad_J4(self):
        loop1 = {'type': 'J3', 'chains': set('')}
        loop2 = {'type': 'J3', 'chains': set('ABCD')}
        assert self.bad_chain_number(loop1) is True
        assert self.bad_chain_number(loop2) is True

    def test_complains_invalid_loop(self):
        loop1 = {'type': 'J4', 'chains': set('ABCD')}
        with pytest.raises(core.InvalidState):
            self.bad_chain_number(loop1)


class RealDataValidationTests(Base):
    loader_class = Loader

    def test_can_detect_if_valid(self):
        loop = self.loop('HL_1GID_001')
        cif = self.loader.cif('1GID')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.status(incomplete, loop) == 1
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'HL_1GID_001',
            'complementary': None,
            'modifications': None,
            'nt_signature': '149, 150, 151, 152, 153, 154',
            'status': 1,
            'loop_release_id': '0.01'
        }

    def test_can_detect_if_has_breaks(self):
        loop = self.loop('HL_1FG0_002')
        cif = self.loader.cif('1FG0')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.has_breaks(loop) is True
        assert self.loader.status(incomplete, loop) == 2
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'HL_1FG0_002',
            'complementary': None,
            'modifications': None,
            'nt_signature': '2136, 2237, 2238, 2239',
            'status': 2,
            'loop_release_id': '0.01'
        }

    def test_can_detect_if_has_modifications(self):
        loop = self.loop('HL_1FCW_015')
        cif = self.loader.cif('1FCW')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.has_modified(loop) is True
        assert self.loader.modified_bases(loop) == ['5MU', 'PSU', '1MA']
        assert self.loader.status(incomplete, loop) == 3
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'HL_1FCW_015',
            'complementary': None,
            'modifications': '5MU, PSU, 1MA',
            'nt_signature': '53, 56, 57, 59, 60, 61',
            'status': 3,
            'loop_release_id': '0.01',
        }

    @pytest.mark.xfail(reason='Bug fix in fr3d prevents this')
    def test_can_detect_if_has_bad_chain_number(self):
        loop = self.loop('HL_1A34_001')
        cif = self.loader.cif('1A34')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.bad_chain_number(loop) is True
        assert self.loader.status(incomplete, loop) == 4
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'HL_1A34_001',
            'complementary': None,
            'modifications': None,
            'nt_signature': '9, 10, 1',
            'status': 4,
            'loop_release_id': '0.01',
        }

    def test_can_detect_incomplete_nts(self):
        loop = self.loop('IL_2HOJ_001')
        cif = self.loader.cif('2HOJ')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.has_incomplete_nucleotides(incomplete, loop) is True
        assert self.loader.status(incomplete, loop) == 5
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'IL_2HOJ_001',
            'complementary': None,
            'modifications': None,
            'nt_signature': '18, 19, 20, 45, 47, 48',
            'status': 5,
            'loop_release_id': '0.01'
        }

    def test_can_detect_is_complementary(self):
        loop = self.loop('IL_1GRZ_007')
        cif = self.loader.cif('1GRZ')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.is_complementary(loop) is True
        assert self.loader.status(incomplete, loop) == 6
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'IL_1GRZ_007',
            'complementary': 'CAG,CUG',
            'modifications': None,
            'nt_signature': '232, 233, 234, 240, 241, 242',
            'status': 6,
            'loop_release_id': '0.01',
        }

    @pytest.mark.skip(reason='Bug fix prevents this')
    def test_can_detect_if_more_than_one_symmetry(self):
        loop = self.loop('HL_1DUH_001')
        cif = self.loader.cif('1DUH')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.too_many_sym_ops(loop) is True
        assert self.loader.status(incomplete, loop) == 7
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'HL_1DUH_001',
            'complementary': 'CAG,CUG',
            'modifications': None,
            'nt_signature': '232, 233, 234, 240, 241, 242',
            'status': 7,
            'loop_release_id': '0.01',
        }

    def test_handles_structures_with_bad_flanking(self):
        loop = self.loop('IL_3J9M_017')
        assert loop['endpoints'] == [
            ('3J9M|1|A|A|2033', '3J9M|1|A|U|2035'),
            ('3J9M|1|A|G|2040', '3J9M|1|A|U|2041'),
        ]

    @pytest.mark.skip(reason="no data yet")
    def test_handles_large_chain_break(self):
        loop = self.loop('HL_4W23_015')
        cif = self.loader.cif('4W23')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.has_breaks(loop) is True
        assert self.loader.status(incomplete, loop) == 2
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'HL_4W23_015',
            'complementary': None,
            'modifications': None,
            'nt_signature': '232, 233, 234, 240, 241, 242',
            'status': 2,
            'loop_release_id': '0.01',
        }


class IncompleteResidueTest(Base):
    loader_class = Loader

    @pytest.mark.xfail(reason="No data yet")
    def test_knows_if_not_incomplete(self):
        loop = self.loop('HL_1GID_003')
        cif = self.loader.cif('1GID')
        incomplete = self.loader.incomplete(cif)
        has_inc = self.loader.has_incomplete_nucleotides(incomplete, loop)
        assert has_inc is False
        assert self.loader.status(incomplete, loop) == 1

    def test_can_detect_if_has_incomplete_nts(self):
        loop = self.loop('IL_2HOM_001')
        cif = self.loader.cif('2HOM')
        incomplete = self.loader.incomplete(cif)
        assert self.loader.has_incomplete_nucleotides(incomplete, loop) is True
        assert self.loader.status(incomplete, loop) == 5
        assert self.loader.quality(incomplete, '0.01', loop) == {
            'loop_id': 'IL_2HOM_001',
            'complementary': None,
            'modifications': None,
            'nt_signature': '18, 19, 20, 45, 47, 48',
            'status': 5,
            'loop_release_id': '0.01',
        }

    def test_it_handles_all_current_examples(self):
        ids = [
            'IL_2HOJ_001',
            'IL_2HOJ_003',
            'IL_2HOK_001',
            'IL_2HOK_003',
            'HL_2HOK_002',
            'IL_2HOL_001',
            'IL_2HOL_003',
            'IL_2HOM_001',
            'IL_2HOM_003',
            'IL_2HOO_001',
            'IL_2HOO_003',
            # 'HL_2HOO_001', Is actually a break
        ]

        for loop_id in ids:
            loop = self.loop(loop_id)
            pdb = loop_id.split('_')[1]
            cif = self.loader.cif(pdb)
            partial = self.loader.incomplete(cif)
            inc = self.loader.has_incomplete_nucleotides(partial, loop)
            assert inc is True, loop_id
            assert self.loader.status(partial, loop) == 5, loop_id
            val = self.loader.quality(partial, '0.01', loop)
            val.pop('complementary')  # Ignoring for this test
            ans = {
                'loop_id': loop_id,
                'modifications': None,
                'nt_signature': val['nt_signature'],
                'status': 5,
                'loop_release_id': '0.01',
            }
            assert val == ans, loop_id


class RealDataTest(Base):
    loader_class = Loader

    def data(self, pdb, status_only=True):
        data = self.loader.data(pdb)
        result = {}
        for entry in data:
            current = entry
            if status_only:
                current = entry['status']
            result[entry['loop_id']] = current
        return result

    def test_can_process_1GID(self):
        assert self.data('1GID') == {
            'HL_1GID_001': 1,
            'HL_1GID_002': 1,
            'HL_1GID_003': 5,
            'HL_1GID_004': 1,
            'HL_1GID_005': 1,
            'HL_1GID_006': 5,
            'IL_1GID_001': 5,
            'IL_1GID_002': 1,
            'IL_1GID_003': 1,
            'IL_1GID_004': 1,
            'IL_1GID_005': 5,
            'IL_1GID_006': 1,
            'IL_1GID_007': 1,
            'IL_1GID_008': 5,
            'IL_1GID_009': 1,
            'IL_1GID_010': 1,
            'IL_1GID_011': 5,
            'IL_1GID_012': 5,
            'IL_1GID_013': 1,
            'IL_1GID_014': 1,
            'J3_1GID_001': 1,
            'J3_1GID_002': 1,
        }

    def test_can_process_2HOJ(self):
        assert self.data('2HOJ') == {
            'IL_2HOJ_001': 5,
            'IL_2HOJ_002': 1,
            'IL_2HOJ_003': 5,
            'HL_2HOJ_001': 2,
            'HL_2HOJ_002': 1
        }

    def test_can_process_124D(self):
        assert self.data('124D') == {}

    def test_can_process_2HOO(self):
        assert self.data('2HOO') == {
            'HL_2HOO_001': 2,
            'HL_2HOO_002': 1,
            'IL_2HOO_001': 5,
            'IL_2HOO_002': 1,
            'IL_2HOO_003': 5,
        }

    def test_can_process_3CPW(self):
        assert self.data('3CPW') == {
            'HL_3CPW_001': 1,
            'HL_3CPW_002': 1,
            'HL_3CPW_003': 1,
            'HL_3CPW_004': 2,
            'HL_3CPW_005': 1,
            'HL_3CPW_006': 1,
            'HL_3CPW_007': 1,
            'HL_3CPW_008': 1,
            'HL_3CPW_009': 1,
            'HL_3CPW_010': 1,
            'HL_3CPW_011': 1,
            'HL_3CPW_012': 1,
            'HL_3CPW_013': 1,
            'HL_3CPW_014': 1,
            'HL_3CPW_015': 1,
            'HL_3CPW_016': 1,
            'HL_3CPW_017': 1,
            'HL_3CPW_018': 1,
            'HL_3CPW_019': 1,
            'HL_3CPW_020': 1,
            'HL_3CPW_021': 2,
            'HL_3CPW_022': 1,
            'HL_3CPW_023': 1,
            'HL_3CPW_024': 1,
            'HL_3CPW_025': 1,
            'HL_3CPW_026': 1,
            'HL_3CPW_027': 2,
            'HL_3CPW_028': 1,
            'HL_3CPW_029': 1,
            'HL_3CPW_030': 1,
            'HL_3CPW_031': 1,
            'HL_3CPW_032': 1,
            'HL_3CPW_033': 1,
            'HL_3CPW_034': 1,
            'HL_3CPW_035': 1,
            'HL_3CPW_036': 1,
            'HL_3CPW_037': 1,
            'HL_3CPW_038': 2,
            'HL_3CPW_039': 1,
            'HL_3CPW_040': 1,
            'HL_3CPW_041': 1,
            'HL_3CPW_042': 1,
            'HL_3CPW_043': 1,
            'HL_3CPW_044': 1,
            'HL_3CPW_045': 1,
            'HL_3CPW_046': 1,
            'HL_3CPW_047': 1,
            'HL_3CPW_048': 1,
            'HL_3CPW_049': 1,
            'HL_3CPW_050': 2,
            'HL_3CPW_051': 1,
            'HL_3CPW_052': 1,
            'HL_3CPW_053': 2,
            'HL_3CPW_054': 1,
            'HL_3CPW_055': 1,
            'HL_3CPW_056': 1,
            'HL_3CPW_057': 2,
            'HL_3CPW_058': 1,
            'HL_3CPW_059': 1,
            'HL_3CPW_060': 1,
            'HL_3CPW_061': 1,
            'HL_3CPW_062': 1,
            'HL_3CPW_063': 1,
            'HL_3CPW_064': 1,
            'HL_3CPW_065': 1,
            'HL_3CPW_066': 1,
            'HL_3CPW_067': 1,
            'HL_3CPW_068': 1,
            'HL_3CPW_069': 1,
            'HL_3CPW_070': 1,
            'HL_3CPW_071': 1,
            'HL_3CPW_072': 1,
            'IL_3CPW_001': 1,
            'IL_3CPW_002': 1,
            'IL_3CPW_003': 1,
            'IL_3CPW_004': 1,
            'IL_3CPW_005': 1,
            'IL_3CPW_006': 1,
            'IL_3CPW_007': 1,
            'IL_3CPW_008': 1,
            'IL_3CPW_009': 1,
            'IL_3CPW_010': 1,
            'IL_3CPW_011': 1,
            'IL_3CPW_012': 1,
            'IL_3CPW_013': 1,
            'IL_3CPW_014': 1,
            'IL_3CPW_015': 1,
            'IL_3CPW_016': 1,
            'IL_3CPW_017': 1,
            'IL_3CPW_018': 1,
            'IL_3CPW_019': 1,
            'IL_3CPW_020': 1,
            'IL_3CPW_021': 1,
            'IL_3CPW_022': 1,
            'IL_3CPW_023': 1,
            'IL_3CPW_024': 1,
            'IL_3CPW_025': 1,
            'IL_3CPW_026': 1,
            'IL_3CPW_027': 1,
            'IL_3CPW_028': 1,
            'IL_3CPW_029': 1,
            'IL_3CPW_030': 1,
            'IL_3CPW_031': 1,
            'IL_3CPW_032': 1,
            'IL_3CPW_033': 1,
            'IL_3CPW_034': 1,
            'IL_3CPW_035': 1,
            'IL_3CPW_036': 1,
            'IL_3CPW_037': 1,
            'IL_3CPW_038': 1,
            'IL_3CPW_039': 1,
            'IL_3CPW_040': 1,
            'IL_3CPW_041': 1,
            'IL_3CPW_042': 1,
            'IL_3CPW_043': 1,
            'IL_3CPW_044': 1,
            'IL_3CPW_045': 1,
            'IL_3CPW_046': 1,
            'IL_3CPW_047': 1,
            'IL_3CPW_048': 1,
            'IL_3CPW_049': 1,
            'IL_3CPW_050': 1,
            'IL_3CPW_051': 1,
            'IL_3CPW_052': 1,
            'IL_3CPW_053': 1,
            'IL_3CPW_054': 1,
            'IL_3CPW_055': 1,
            'IL_3CPW_056': 1,
            'IL_3CPW_057': 1,
            'IL_3CPW_058': 1,
            'IL_3CPW_059': 1,
            'IL_3CPW_060': 1,
            'IL_3CPW_061': 1,
            'IL_3CPW_062': 1,
            'IL_3CPW_063': 1,
            'IL_3CPW_064': 1,
            'IL_3CPW_065': 1,
            'IL_3CPW_066': 1,
            'IL_3CPW_067': 1,
            'IL_3CPW_068': 1,
            'IL_3CPW_069': 1,
            'IL_3CPW_070': 1,
            'IL_3CPW_071': 1,
            'IL_3CPW_072': 1,
            'IL_3CPW_073': 1,
            'IL_3CPW_074': 1,
            'IL_3CPW_075': 1,
            'IL_3CPW_076': 1,
            'IL_3CPW_077': 1,
            'IL_3CPW_078': 1,
            'IL_3CPW_079': 1,
            'IL_3CPW_080': 1,
            'IL_3CPW_081': 1,
            'IL_3CPW_082': 1,
            'IL_3CPW_083': 1,
            'IL_3CPW_084': 1,
            'IL_3CPW_085': 1,
            'IL_3CPW_086': 1,
            'IL_3CPW_087': 1,
            'IL_3CPW_088': 1,
            'IL_3CPW_089': 1,
            'IL_3CPW_090': 1,
            'IL_3CPW_091': 1,
            'IL_3CPW_092': 1,
            'IL_3CPW_093': 1,
            'IL_3CPW_094': 1,
            'IL_3CPW_095': 1,
            'IL_3CPW_096': 1,
            'IL_3CPW_097': 1,
            'IL_3CPW_098': 1,
            'IL_3CPW_099': 1,
            'IL_3CPW_100': 1,
            'IL_3CPW_101': 1,
            'IL_3CPW_102': 1,
            'IL_3CPW_103': 1,
            'IL_3CPW_104': 1,
            'IL_3CPW_105': 1,
            'J3_3CPW_001': 1,
            'J3_3CPW_002': 1,
            'J3_3CPW_003': 1,
            'J3_3CPW_004': 1,
            'J3_3CPW_005': 1,
            'J3_3CPW_006': 1,
            'J3_3CPW_007': 1,
            'J3_3CPW_008': 1,
            'J3_3CPW_009': 1,
            'J3_3CPW_010': 1,
        }

    def test_can_process_2HOM(self):
        assert self.data('2HOM') == {
            'HL_2HOM_001': 2,
            'HL_2HOM_002': 1,
            'IL_2HOM_001': 5,
            'IL_2HOM_002': 1,
            'IL_2HOM_003': 5,
        }

    def test_can_process_1S72(self):
        assert self.data('1S72') == {
            'HL_1S72_001': 1,
            'HL_1S72_002': 1,
            'HL_1S72_003': 1,
            'HL_1S72_004': 2,
            'HL_1S72_005': 1,
            'HL_1S72_006': 1,
            'HL_1S72_007': 1,
            'HL_1S72_008': 1,
            'HL_1S72_009': 1,
            'HL_1S72_010': 1,
            'HL_1S72_011': 1,
            'HL_1S72_012': 1,
            'HL_1S72_013': 1,
            'HL_1S72_014': 1,
            'HL_1S72_015': 1,
            'HL_1S72_016': 1,
            'HL_1S72_017': 1,
            'HL_1S72_018': 3,
            'HL_1S72_019': 1,
            'HL_1S72_020': 1,
            'HL_1S72_021': 2,
            'HL_1S72_022': 1,
            'HL_1S72_023': 1,
            'HL_1S72_024': 1,
            'HL_1S72_025': 1,
            'HL_1S72_026': 1,
            'HL_1S72_027': 2,
            'HL_1S72_028': 1,
            'HL_1S72_029': 1,
            'HL_1S72_030': 1,
            'HL_1S72_031': 1,
            'HL_1S72_032': 1,
            'HL_1S72_033': 1,
            'HL_1S72_034': 1,
            'HL_1S72_035': 1,
            'HL_1S72_036': 1,
            'HL_1S72_037': 1,
            'HL_1S72_038': 2,
            'HL_1S72_039': 1,
            'HL_1S72_040': 1,
            'HL_1S72_041': 1,
            'HL_1S72_042': 1,
            'HL_1S72_043': 1,
            'HL_1S72_044': 1,
            'HL_1S72_045': 1,
            'HL_1S72_046': 1,
            'HL_1S72_047': 1,
            'HL_1S72_048': 1,
            'HL_1S72_049': 1,
            'HL_1S72_050': 2,
            'HL_1S72_051': 1,
            'HL_1S72_052': 1,
            'HL_1S72_053': 2,
            'HL_1S72_054': 1,
            'HL_1S72_055': 1,
            'HL_1S72_056': 1,
            'HL_1S72_057': 2,
            'HL_1S72_058': 1,
            'HL_1S72_059': 1,
            'HL_1S72_060': 1,
            'HL_1S72_061': 1,
            'HL_1S72_062': 1,
            'HL_1S72_063': 1,
            'HL_1S72_064': 3,
            'HL_1S72_065': 1,
            'HL_1S72_066': 1,
            'HL_1S72_067': 1,
            'HL_1S72_068': 1,
            'HL_1S72_069': 1,
            'HL_1S72_070': 1,
            'HL_1S72_071': 1,
            'HL_1S72_072': 1,
            'IL_1S72_001': 1,
            'IL_1S72_002': 1,
            'IL_1S72_003': 1,
            'IL_1S72_004': 1,
            'IL_1S72_005': 1,
            'IL_1S72_006': 1,
            'IL_1S72_007': 1,
            'IL_1S72_008': 1,
            'IL_1S72_009': 1,
            'IL_1S72_010': 1,
            'IL_1S72_011': 1,
            'IL_1S72_012': 1,
            'IL_1S72_013': 1,
            'IL_1S72_014': 1,
            'IL_1S72_015': 1,
            'IL_1S72_016': 1,
            'IL_1S72_017': 1,
            'IL_1S72_018': 1,
            'IL_1S72_019': 1,
            'IL_1S72_020': 1,
            'IL_1S72_021': 1,
            'IL_1S72_022': 1,
            'IL_1S72_023': 1,
            'IL_1S72_024': 1,
            'IL_1S72_025': 1,
            'IL_1S72_026': 1,
            'IL_1S72_027': 1,
            'IL_1S72_028': 1,
            'IL_1S72_029': 1,
            'IL_1S72_030': 1,
            'IL_1S72_031': 1,
            'IL_1S72_032': 1,
            'IL_1S72_033': 1,
            'IL_1S72_034': 1,
            'IL_1S72_035': 1,
            'IL_1S72_036': 1,
            'IL_1S72_037': 1,
            'IL_1S72_038': 1,
            'IL_1S72_039': 1,
            'IL_1S72_040': 1,
            'IL_1S72_041': 1,
            'IL_1S72_042': 1,
            'IL_1S72_043': 1,
            'IL_1S72_044': 1,
            'IL_1S72_045': 1,
            'IL_1S72_046': 1,
            'IL_1S72_047': 1,
            'IL_1S72_048': 1,
            'IL_1S72_049': 1,
            'IL_1S72_050': 1,
            'IL_1S72_051': 1,
            'IL_1S72_052': 1,
            'IL_1S72_053': 1,
            'IL_1S72_054': 1,
            'IL_1S72_055': 1,
            'IL_1S72_056': 1,
            'IL_1S72_057': 1,
            'IL_1S72_058': 1,
            'IL_1S72_059': 1,
            'IL_1S72_060': 1,
            'IL_1S72_061': 1,
            'IL_1S72_062': 1,
            'IL_1S72_063': 1,
            'IL_1S72_064': 1,
            'IL_1S72_065': 1,
            'IL_1S72_066': 1,
            'IL_1S72_067': 1,
            'IL_1S72_068': 1,
            'IL_1S72_069': 1,
            'IL_1S72_070': 1,
            'IL_1S72_071': 1,
            'IL_1S72_072': 1,
            'IL_1S72_073': 1,
            'IL_1S72_074': 1,
            'IL_1S72_075': 1,
            'IL_1S72_076': 1,
            'IL_1S72_077': 1,
            'IL_1S72_078': 1,
            'IL_1S72_079': 1,
            'IL_1S72_080': 1,
            'IL_1S72_081': 1,
            'IL_1S72_082': 1,
            'IL_1S72_083': 1,
            'IL_1S72_084': 1,
            'IL_1S72_085': 1,
            'IL_1S72_086': 1,
            'IL_1S72_087': 1,
            'IL_1S72_088': 1,
            'IL_1S72_089': 1,
            'IL_1S72_090': 1,
            'IL_1S72_091': 1,
            'IL_1S72_092': 1,
            'IL_1S72_093': 1,
            'IL_1S72_094': 1,
            'IL_1S72_095': 1,
            'IL_1S72_096': 1,
            'IL_1S72_097': 1,
            'IL_1S72_098': 1,
            'IL_1S72_099': 1,
            'IL_1S72_100': 1,
            'IL_1S72_101': 1,
            'IL_1S72_102': 1,
            'IL_1S72_103': 1,
            'IL_1S72_104': 1,
            'J3_1S72_001': 1,
            'J3_1S72_002': 1,
            'J3_1S72_003': 1,
            'J3_1S72_004': 1,
            'J3_1S72_005': 1,
            'J3_1S72_006': 1,
            'J3_1S72_007': 1,
            'J3_1S72_008': 1,
            'J3_1S72_009': 1,
            'J3_1S72_010': 1,
        }
