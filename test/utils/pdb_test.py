from datetime import date
from unittest import TestCase

from pymotifs.utils.pdb import RnaPdbsHelper


class RnaPdbsHelperTest(TestCase):

    def setUp(self):
        self.helper = RnaPdbsHelper()

    def test_can_find_all_pdbs(self):
        val = self.helper()
        self.assertTrue(len(val) >= 3000)

    def test_can_find_all_before_date(self):
        val = self.helper(dates=(None, date(2014, 12, 11)))
        self.assertEquals(2712, len(val))

    def test_can_find_all_after_date(self):
        val = self.helper(dates=(date(2014, 12, 11), None))
        self.assertTrue(len(val) > 300)

    def test_can_find_all_between_dates(self):
        val = self.helper(dates=(date(2014, 12, 11), date(2014, 12, 20)))
        self.assertEquals(['4WAL', '4WAN'], val)

    def test_can_find_between_dates_with_nothing(self):
        val = self.helper(dates=(date(2014, 12, 11), date(2014, 12, 16)))
        self.assertEquals([], val)

    def test_fails_given_invalid_dates(self):
        self.assertRaises(Exception, self.helper,
                          dates=('bob', date(2014, 12, 15)))
