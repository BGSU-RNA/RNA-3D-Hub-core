from test import StageTest
from nose import SkipTest

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

    def test_process_gets_all_at_once(self):
        val = ['A', 'B', 'C', 'D']
        self.loader.process = lambda a: self.assertEquals(tuple(val), a)
        self.loader.data(val)

    def test_to_process_creates_a_tuple(self):
        val = ['A', 'C', 'c', 'D']
        self.assertEquals(tuple(['A', 'C', 'C', 'D']),
                          self.loader.to_process(val))

    def test_data_gets_all_at_once(self):
        raise SkipTest()
        val = ['A', 'B', 'C', 'D']

        def data(args):
            self.assertEquals(tuple(val), args)
        self.loader.data = data
        self.loader.data(val)

    def test_mark_processed_gets_all_at_once(self):
        val = [1, 3]
        self.loader.mark_processed = lambda a: self.assertEquals(val, a)
        self.loader.mark_processed(val)


class MarkingProcessedTest(StageTest):
    loader_class = Mass

    def test_marks_each_done_individually(self):
        raise SkipTest()
