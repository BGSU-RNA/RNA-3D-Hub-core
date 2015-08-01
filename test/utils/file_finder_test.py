import os
import unittest

from pymotifs.utils import CifFileFinder
from pymotifs.utils import MissingFileException


class CifFileFinderTest(unittest.TestCase):
    def setUp(self):
        path = os.path.join(os.path.normpath("."), 'FR3D')
        self.cif = CifFileFinder({'locations': {'fr3d_root': path}})

    def test_can_find_a_cif_file(self):
        val = self.cif("1GID")
        ans = os.path.abspath("./FR3D/PDBFiles/1GID.cif")
        self.assertEqual(ans, val)

    def test_fails_for_missing_file(self):
        self.assertRaises(MissingFileException, self.cif, "bob")
