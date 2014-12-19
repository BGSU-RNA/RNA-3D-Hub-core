import unittest as ut

from pymotifs.pdbs import info
from pymotifs.models import Session


class InfoTest(ut.TestCase):
    def setUp(self):
        self.loader = info.PdbInfoLoader({}, Session)

    def test_knows_if_data_is_present(self):
        self.assertTrue(self.loader.has_data('2AW7'))

    def test_knows_if_data_is_missing(self):
        self.assertFalse(self.loader.has_data('bob'))

    # def test_it_can_remove_data(self):
    #     self.loader.remove('2AW7')
    #     self.assertFalse(self.loader.has_data('2AW7'))

    def test_it_can_create_data(self):
        pass
