"""



"""

import unittest


import UnitIdLoader as loader
import models


class TestUnitIdLoader(unittest.TestCase):

    def setUp(self):
        self.pdbs = ['2Z4L'] # need to choose a file with both au and ba1
        L = loader.UnitIdLoader()
        L.import_unit_ids(self.pdbs, recalculate=True)

    def test_import(self):
        ids = models.session.query(models.PdbUnitIdCorrespondence).\
                             filter_by(pdb='2Z4L').\
                             all()
        self.assertEqual(len(ids), 14256)

    def tearDown(self):
        models.session.query(models.PdbUnitIdCorrespondence).delete()


if __name__ == '__main__':
    unittest.main()