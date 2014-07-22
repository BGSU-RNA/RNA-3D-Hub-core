import os
import unittest

from chain_breaks import ChainBreakFinder


class ChainBreakFinderTest(unittest.TestCase):

    def setUp(self):
        self.finder = ChainBreakFinder()

    def test_finds_chain_breaks(self):
        breaks = self.finder(os.path.normpath("../FR3D/PDBFiles/2AW7.cif"))
        self.assertEqual(22, len(breaks))
