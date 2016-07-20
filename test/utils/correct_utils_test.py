import pytest

from test import StageTest

from pymotifs import core
from pymotifs.utils.correct_units import Unit
from pymotifs.utils.correct_units import Correcter


class NormalizedMappingTest(StageTest):
    loader_class = Correcter

    def test_can_create_complete_mapping(self):
        mapping = self.loader.normalized_mapping('157D')
        assert len(mapping) == 121

    def test_can_create_a_normalized_mapping(self):
        mapping = self.loader.normalized_mapping('157D')
        assert self.loader.as_unit('157D|1|A|G|10') in mapping

    def test_creates_mapping_for_sym_ops(self):
        mapping = self.loader.normalized_mapping('1A34')
        assert len(mapping) == 698
        assert self.loader.as_unit('1A34|1|A|ALA|25||||P_1') in mapping



class UtilsTest(StageTest):
    loader_class = Correcter

    def test_can_turn_a_unit_id_to_unit(self):
        unit = Unit(model=2, chain='A', component_number=6,
                    insertion_code=None, alt_id=None, symmetry='6_765')
        assert self.loader.as_unit('1Z7F|2|A|C|6||||6_765') == unit


class CorrectorsTest(StageTest):
    loader_class = Correcter

    def test_will_set_model_to_one_if_needed(self):
        unit = Unit(model=2, chain='A', component_number=6,
                    insertion_code=None, alt_id=None, symmetry='6_765')
        ans = Unit(model=1, chain='A', component_number=6, insertion_code=None,
                   alt_id=None, symmetry='6_765')
        assert self.loader.correct_model(unit) == ans

    def test_will_set_model_to_one_if_not_needed(self):
        unit = Unit(model=2, chain='A', component_number=6,
                    insertion_code=None, alt_id=None, symmetry='1_555')
        assert self.loader.correct_model(unit) is None

    def test_will_correct_alt_id_if_needed(self):
        unit = Unit(model=2, chain='A', component_number=6,
                    insertion_code=None, alt_id=None, symmetry='6_765')
        ans = Unit(model=2, chain='A', component_number=6,
                   insertion_code=None, alt_id='a', symmetry='6_765')
        assert self.loader.correct_alt_id(unit) == ans

    def will_not_correct_alt_id_if_not_needed(self):
        unit = Unit(model=2, chain='A', component_number=6,
                    insertion_code=None, alt_id='b', symmetry='6_765')
        assert self.loader.correct_alt_id(unit) is None

    def test_will_correct_alt_id_and_model_if_needed(self):
        unit = Unit(model=2, chain='A', component_number=6,
                    insertion_code=None, alt_id=None, symmetry='6_765')
        ans = Unit(model=1, chain='A', component_number=6,
                   insertion_code=None, alt_id='a', symmetry='6_765')
        assert self.loader.correct_model_and_alt_id(unit) == ans

    def test_will_not_correct_alt_id_and_model_if_not_needed(self):
        unit = Unit(model=1, chain='A', component_number=6,
                    insertion_code=None, alt_id='b', symmetry='6_765')
        assert self.loader.correct_model_and_alt_id(unit) is None

    def test_can_do_nothing_to_unit(self):
        unit = Unit(model=2, chain='A', component_number=6, insertion_code=None,
                    alt_id=None, symmetry='6_765')
        assert self.loader.correct_nothing(unit) == unit


class CorrectingTest(StageTest):
    loader_class = Correcter

    def build(self, *unit_ids):
        mapping = {}
        for unit_id in unit_ids:
            mapping[self.loader.as_unit(unit_id)] = unit_id
        return mapping

    def correct(self, unit_id):
        mapping = self.build('1Z7F|1|A|A|10', '1Z7F|1|A|A|10||||6_765',
                             '1Z7F|1|A|A|11', '1Z7F|1|A|A|11||||6_765',
                             '1Z7F|1|A|A|12||a', '1Z7F|1|A|A|12||a||6_765',
                             )
        unit = self.loader.as_unit(unit_id)
        corrected = self.loader.correct(mapping, unit)
        if corrected:
            return mapping[corrected]
        return None

    def test_it_will_not_correct_if_unneeded(self):
        assert self.correct('1Z7F|1|A|A|10') == '1Z7F|1|A|A|10'

    def test_it_will_correct_model(self):
        assert self.correct('1Z7F|2|A|A|11||||6_765') == '1Z7F|1|A|A|11||||6_765'

    def test_it_will_correct_alt_id(self):
        assert self.correct('1Z7F|1|A|A|12') == '1Z7F|1|A|A|12||a'

    def test_it_will_correct_model_and_alt_id(self):
        assert self.correct('1Z7F|2|A|A|12||||6_765') == '1Z7F|1|A|A|12||a||6_765'

    def test_will_fail_if_cannot_correct(self):
        self.correct('bob') == None


class RealDataTest(StageTest):
    loader_class = Correcter

    @pytest.mark.skip()
    def test_1Z7F(self):
        corrected = self.loader(['1Z7F|2|A|A|5||||6_765',
                                 '1Z7F|2|A|C|6||||6_765',
                                 '1Z7F|2|A|U|7||||6_765',
                                 '1Z7F|1|A|A|10',
                                 '1Z7F|1|A|A|11',
                                 '1Z7F|1|A|U|12'])
        assert corrected == ['1Z7F|1|A|A|5||||6_765',
                             '1Z7F|1|A|C|6||||6_765',
                             '1Z7F|1|A|U|7||||6_765',
                             '1Z7F|1|A|A|10',
                             '1Z7F|1|A|A|11',
                             '1Z7F|1|A|U|12']
