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


class ModifiedBasesTest(StageTest):
    loader_class = Loader

    def test_knows_if_unmodified_bases(self):
        valid = {'units': ['A', 'C', 'G', 'U']}
        assert [] == self.loader.modified_bases(valid)
        assert False is self.loader.has_modified(valid)

    def test_knows_if_has_modified_bases(self):
        invalid = {'units': ['A', 'C', 'G', 'I']}
        assert ['I'] == self.loader.modified_bases(invalid)
        assert True is self.loader.has_modified(invalid)


class ComplementarySequenceTest(StageTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_can_detect_if_has_non_cWW_pairs(self):
        pass

    @pytest.mark.skip()
    def test_can_detect_if_is_complementary(self):
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

    @pytest.mark.skip()
    def test_can_detect_if_has_breaks(self):
        loop = self.loop('HL_1FG0_002')
        assert self.loader.has_breaks(loop) is True
        assert self.loader.status(loop) == 2

    @pytest.mark.skip()
    def test_can_detect_if_has_modifications(self):
        loop = self.loop('HL_1FCW_015')
        assert self.loader.has_modified(loop) == '5MU, PSU, 1MA'
        assert self.loader.status(loop) == 3

    @pytest.mark.skip()
    def test_can_detect_if_has_bad_chain_number(self):
        loop = self.loop('HL_1FEU_001')
        assert self.loader.has_bad_chain_number(loop) is True
        assert self.loader.status(loop) == 4

    @pytest.mark.skip()
    def test_can_detect_if_has_incomplete_nts(self):
        loop = self.loop('IL_2HOM_001')
        assert self.loader.has_incomplete_nucleotides(loop) is True
        assert self.loader.status(loop) == 5

    @pytest.mark.skip()
    def test_can_detect_if_is_complemenatry(self):
        loop = self.loop('IL_1GRZ_007')
        assert self.loader.is_complementary(loop) == 'CAG,CUG'
        assert self.loader.status(loop) == 6
