from pymotifs import core
import pymotifs.quality.utils as ut

import pytest

from test import StageTest

from fr3d.unit_ids import decode


class ParserTest(StageTest):
    loader_class = ut.Utils
    filename = None

    @classmethod
    def setUpClass(cls):
        super(ParserTest, cls).setUpClass()
        with open(cls.filename, 'rb') as raw:
            cls.parser = ut.Parser(raw.read())

    def setUp(self):
        super(ParserTest, self).setUp()
        self.parser = self.__class__.parser

    def mapping(self, pdb):
        return self.loader.unit_mapping(pdb)


class CoreRsrParserTest(ParserTest):
    filename = 'test/files/validation/4v7w_validation.xml.gz'

    def test_can_get_tree_from_gz_content(self):
        self.assertTrue(self.parser.root)

    def test_can_detect_has_rsr(self):
        self.assertTrue(self.parser.has_rsr())

    def test_can_detect_has_dcc(self):
        assert self.parser.has_dcc() is True

    def test_can_map_all_nts(self):
        mapping = self.mapping('4V7W')
        assert len(list(self.parser.nts(mapping))) == 20917

    def test_fails_if_cannot_map_all(self):
        mapping = self.mapping('4V7W')
        mapping.pop(mapping.keys()[0])
        with pytest.raises(core.InvalidState):
            list(self.parser.nts(mapping))

    def test_can_generate_nt_level_data(self):
        mapping = self.mapping('4V7W')
        assert list(self.parser.nts(mapping))[0] == {
            'id': '4V7W|1|AA|U|5',
            'real_space_r': 0.218,
            'real_space_r_z_score': 0.26,
        }

    def test_can_generate_structure_level_data(self):
        assert self.parser.entity() == {
            'pdb_id': '4V7W',
            'percent_rsrz_outliers': 6.76,
            'absolute_percentile_percent_rsrz_outliers': 17.0,
            'relative_percentile_percent_rsrz_outliers': 4.3,
            'clashscore': 40.92,
            'relative_percentile_clashscore': 16.5,
            'absolute_percentile_clashscore': 3.1,
            'percent_rota_outliers': 17.72,
            'absolute_percentile_percent_rota_outliers': 2.5,
            'relative_percentile_percent_rota_outliers': 13.8,
            'md5': 'ad9cd539ce3e7f8c83d1fa706bf3c79a',
        }


class MissingRsRParserTest(ParserTest):
    filename = 'test/files/validation/1j5e_validation.xml.gz'

    def test_can_tell_has_no_rsr(self):
        assert self.parser.has_rsr() is False

    def test_can_tell_has_no_dcc(self):
        assert self.parser.has_dcc() is False


class MissingDataTest(ParserTest):
    filename = 'test/files/validation/1fjg_validation.xml.gz'

    def setUp(self):
        super(MissingDataTest, self).setUp()
        self.parser = self.__class__.parser
        mapping = self.mapping('1FJG')
        self.nts = list(self.parser.nts(mapping))

    def chain_of(self, uid):
        return decode(uid)['chain']

    def test_parses_all_rna_data(self):
        val = [nt for nt in self.nts if self.chain_of(nt['id']) == 'A']
        self.assertEquals(len(val), 1603)


class AltIdParsingTest(ParserTest):
    filename = 'test/files/validation/1vy4_validation.xml.gz'

    def setUp(self):
        super(AltIdParsingTest, self).setUp()
        self.parser = self.__class__.parser
        mapping = self.mapping('1VY4')
        self.nts = list(self.parser.nts(mapping))

    def alt_of(self, uid):
        return decode(uid)['chain']

    def test_can_generate_ids_using_alt_ids(self):
        val = [nt for nt in self.nts if self.alt_of(nt['id'])]
        self.assertTrue(2, len(val))


class UnusualUnitsTest(ParserTest):
    filename = 'test/files/validation/2uua_validation.xml.gz'

    def setUp(self):
        super(UnusualUnitsTest, self).setUp()
        self.parser = self.__class__.parser
        mapping = self.mapping('2UUA')
        self.nts = list(self.parser.nts(mapping))

    def test_can_generate_unit_ids_for_all_units(self):
        assert len(self.nts) == 4129

    @pytest.mark.skip()
    def test_can_generate_unit_id_for_hard_unit(self):
        data = {
            'real_space_r': 0.32,
            'id': {
                'component_id': 'PAR',
                'chain': 'Z',
                'insertion_code': None,
                'component_number': 1,
                'alt_id': None,
                'model': 1,
                'pdb': '2UUA'
            }
        }


class UtilsTest(StageTest):
    loader_class = ut.Utils

    def test_can_get_filename(self):
        ans = '/Users/bsweene/dotfiles/personal/leontis/hub-core/MotifAtlas/quality/validation-reports/1FJG.xml.gz'
        assert self.loader.filename('1FJG') == ans

    def test_knows_if_has_no_data(self):
        assert self.loader.has_no_data('0FJG') is True

    @pytest.mark.skip()
    def test_knows_has_no_data_if_empty(self):
        assert self.loader.has_no_data('') is True

    @pytest.mark.skip()
    def test_can_list_known_reports(self):
        assert self.loader.known() == []

    @pytest.mark.skip()
    def test_can_list_reports_with_data(self):
        assert self.loader.known(has_data=True) == []

    @pytest.mark.skip()
    def test_can_list_reports_without_data(self):
        assert self.loader.known(has_data=False) == []

    def test_can_create_a_unit_mapping(self):
        val = self.loader.unit_mapping('124D')
        assert len(val) == 16
        assert val[('A', 4, None, None)] == ['124D|1|A|DA|4']

    def test_can_create_mapping_with_sym_ops(self):
        val = self.loader.unit_mapping('1A34')
        ans = ['1A34|1|C|U|7||||P_1', '1A34|1|C|U|7||||P_P']
        assert len(val) == 349
        assert val[('C', 7, None, None)] == ans

    def test_can_create_mapping_with_alt_ids(self):
        val = self.loader.unit_mapping('1A34')
        ans = ['1A34|1|A|CYS|157||A||P_1', '1A34|1|A|CYS|157||A||P_P']
        assert len(val) == 349
        assert val[('A', 157, None, 'A')] == ans

    def test_can_create_mapping_with_insertion_codes(self):
        val = self.loader.unit_mapping('1FJG')
        ans = ['1FJG|1|A|A|1030|||D']
        assert len(val) == 4018
        assert val[('A', 1030, 'D', None)] == ans
