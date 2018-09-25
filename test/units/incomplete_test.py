from test import StageTest

from pymotifs.units.incomplete import Loader
from pymotifs.units.incomplete import Entry


class QueryTest(StageTest):
    loader_class = Loader

    def test_can_determine_if_has_data(self):
        assert self.loader.has_data('1FJG') is True

    def test_can_determine_if_has_no_data(self):
        assert self.loader.has_data('0FJG') is False


class MissingTest(StageTest):
    loader_class = Loader

    def test_can_find_all_missing_or_incomplete(self):
        val = self.loader.missing_keys('1FJG')
        assert len(val) == 154

    def test_includes_incomplete_entries(self):
        val = self.loader.missing_keys('1FJG')
        assert Entry('1FJG', 1, 'A', 1534, 'A', None, None) in val
        assert Entry('1FJG', 1, 'A', 5, 'U', None, None) in val
        assert Entry('1FJG', 1, 'J', 100, 'THR', None, None) in val
        assert len(val) == 154
