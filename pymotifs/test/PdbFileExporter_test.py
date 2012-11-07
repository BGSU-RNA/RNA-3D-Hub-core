"""



"""

import unittest
import os


import PdbFileExporter as exporter
import models


class TestPdbFileExporter(unittest.TestCase):

    def setUp(self):
        E = exporter.PdbFileExporter()
        self.output = os.path.join(os.getcwd(), 'TestExportInteractions.csv.gz')
        E.export_interactions(self.output, pdb_ids=['1J5E'])

    def test_export(self):
        """
            Test for the presence of the file. Under test environment, the file
            will contain only header. More thorough testing requires setting up
            pdb_coordinates and other tables.
        """
        self.assertTrue( os.path.exists(self.output) )

    def tearDown(self):
        os.remove(self.output)


if __name__ == '__main__':
    unittest.main()