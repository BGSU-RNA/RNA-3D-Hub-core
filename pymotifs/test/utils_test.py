import os
import unittest

from utils import CifFileFinder
from utils import MissingFileException


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
