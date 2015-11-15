import os
import unittest

from test import CONFIG

from pymotifs.utils import CifFileFinder
from pymotifs.utils import MissingFileException


class CifFileFinderTest(unittest.TestCase):
    def setUp(self):
        self.cif = CifFileFinder(CONFIG)

    def test_can_find_a_cif_file(self):
        val = self.cif("1GID")
        ans = os.path.abspath("./FR3D/PDBFiles/1GID.cif")
        self.assertEqual(ans, val)

    def test_fails_for_missing_file(self):
        self.assertRaises(MissingFileException, self.cif, "bob")
