from test import StageTest

from pymotifs.units.incomplete import Loader
from pymotifs.units.incomplete import Entry


class QueryTest(StageTest):
    loader_class = Loader

    def test_can_determine_if_has_data(self):
        assert self.loader.has_data('1FJG') is True

    def test_can_determine_if_has_no_data(self):
        assert self.loader.has_data('0FJG') is False


class MappingTest(StageTest):
    loader_class = Loader

    def test_can_build_complete_mapping(self):
        val = self.loader.mapping('124D')
        assert len(val) == 16
        key = Entry(1, 'A', 4, 'DA', None, None)
        assert val[key] == set(['124D|1|A|DA|4'])

    def test_builds_mapping_without_symmetries(self):
        val = self.loader.mapping('1A34')
        assert len(val) == 349
        # 1A34|1|A|ALA|137||||P_P
        key = Entry(1, 'A', 137, 'ALA', None, None)
        assert val[key] == set(['1A34|1|A|ALA|137||||P_P',
                                '1A34|1|A|ALA|137||||P_1'])


class MissingTest(StageTest):
    loader_class = Loader

    def test_can_find_all_missing_or_incomplete(self):
        val = self.loader.missing_keys('1FJG')
        assert len(val) == 154

    def test_includes_incomplete_entries(self):
        val = self.loader.missing_keys('1FJG')
        assert Entry(1, 'A', 1534, 'A', None, None) in val
        assert Entry(1, 'A', 5, 'U', None, None) in val
        assert Entry(1, 'J', 100, 'THR', None, None) in val


class DataTest(StageTest):
    loader_class = Loader

    def test_can_create_entries_for_all_incomplete(self):
        data = self.loader.data('1FJG')
        val = sorted(d.unit_id for d in data)
        assert val == ['1FJG|1|A|A|1534', '1FJG|1|A|U|5', '1FJG|1|J|THR|100']
