from test import StageTest

from pymotifs.core import MassLoader


class Mass(MassLoader):
    def data(self, *args, **kwargs):
        pass


class ToProecssTest(StageTest):
    loader_class = Mass

    def test_creates_a_tuple_of_all(self):
        self.assertEquals(('A', 'B'), self.loader.to_process(['A', 'B']))

    def test_up_cases_all_entries(self):
        self.assertEquals(('A', 'B'), self.loader.to_process(['a', 'b']))


class BasicTests(StageTest):
    loader_class = Mass

    def tests_returns_false_for_has_data(self):
        self.assertFalse(self.loader.has_data('bob'))


class GeneratingDataTest(StageTest):
    loader_class = Mass

    def test_processes_all_at_once(self):
        val = ['A', 'B', 'C', 'D']
        self.loader.data = lambda args: self.assertEquals(tuple(val), args)
        self.loader.data(val)
