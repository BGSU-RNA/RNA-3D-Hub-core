import pytest

from pymotifs import core

from test import StageTest

from pymotifs.loops.quality import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_can_detect_if_has_data(self):
        self.assertTrue(self.loader.has_data('1GID'))

    def test_can_detect_if_has_no_data(self):
        self.assertFalse(self.loader.has_data('0GID'))


class ComplementarySequenceTest(StageTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_can_detect_if_has_non_cWW_pairs(self):
        pass

    @pytest.mark.skip()
    def test_can_find_if_has_cWW_pairs(self):
        pass

    def test_can_find_if_complemenatry_sequence(self):
        valid = {'units': ['A', 'C', 'C', 'G']}
        invalid1 = {'units': ['A', 'C', 'G', 'U']}
        invalid2 = {'units': ['A', 'C', 'A', 'G', 'U']}
        assert self.loader.complementary_sequence(valid) is False
        assert self.loader.complementary_sequence(invalid1) is True
        assert self.loader.complementary_sequence(invalid2) is True


class ChainNumberValidationTests(StageTest):
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


class RealDataValidationTests(StageTest):
    loader_class = Loader

    def loop(self, loop_id):
        pdb = loop_id.split('_')[1]
        for loop in self.loader.loops(pdb):
            if loop['id'] == loop_id:
                return loop
        self.fail("Could not find loop")

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

    @pytest.mark.skip()
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
            'nt_signature': '23, 24, 25, 29, 30, 31, 32, 33, 34, 35, 36',
            'status': 5,
            'loop_release_id': '0.01',
        }

    def test_can_detect_if_is_complemenatry(self):
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


class IncompleteResidueTest(StageTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_it_handles_all_current_examples(self):
        ids = [
            'HL_1DUH_001',
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
            'HL_2HOO_001',
            'HL_1H3E_001',
        ]

        for loop_id in ids:
            loop = self.loop(loop_id)
            pdb = loop_id.split('_')[1]
            cif = self.loader.cif(pdb)
            partial = self.loader.incomplete(cif)
            inc = self.loader.has_incomplete_nucleotides(partial, loop)
            assert inc is True
            assert self.loader.status(partial, loop) == 5
            val = self.loader.quality(partial, '0.01', loop)
            assert val == {
                'loop_id': loop_id,
                'complementary': None,
                'modifications': None,
                'nt_signature': val['nt_signature'],
                'status': 5,
                'loop_release_id': '0.01',
            }
