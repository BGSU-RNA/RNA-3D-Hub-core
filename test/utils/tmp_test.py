import unittest as ut

from pymotifs.utils import tmp


class StoringTest(ut.TestCase):
    def test_can_store_and_load_data(self):
        tmp.store('test', [1])
        self.assertEqual([1], tmp.load('test'))

    def test_gives_none_for_unknown_group(self):
        self.assertEqual(None, tmp.load('bob'))

    def test_cleans_unknown_group(self):
        self.assertTrue(tmp.cleanup('bob'))

    def test_can_cleanup(self):
        tmp.store('test', [1])
        tmp.cleanup('test')
        self.assertEqual(None, tmp.load('test'))
