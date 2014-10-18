import os
import unittest

from utils import CifFileFinder
from utils import MissingFileException
import utils


class CifFileFinderTest(unittest.TestCase):
    def setUp(self):
        path = os.path.normpath("..")
        self.cif = CifFileFinder({'locations': {'fr3d_root': path}})

    def test_can_find_a_cif_file(self):
        val = self.cif("2AW7")
        ans = os.path.normpath("../FR3D/PDBFiles/2AW7.cif")
        self.assertEqual(ans, val)

    def test_fails_for_missing_file(self):
        self.assertRaises(MissingFileException, self.cif, "bob")


class FailingRetry(utils.RetryHelper):
    def action(*args, **kwargs):
        raise ValueError("Expected")


class RetryHelperTest(unittest.TestCase):
    def test_can_allow_failure(self):
        val = FailingRetry(allow_fail=True)()
        self.assertTrue(val is None)

    def test_will_fail_with_failure(self):
        self.assertRaises(utils.RetryFailedException, FailingRetry())
