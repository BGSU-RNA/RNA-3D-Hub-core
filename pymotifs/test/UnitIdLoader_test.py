"""



"""

import unittest


import UnitIdLoader as loader
import models


class TestUnitIdLoader(unittest.TestCase):

    def setUp(self):
        # 2A34 has both au and ba1, 1A34 is a viral icosahedral structure
        self.pdbs = ['2Z4L', '1A34']
        L = loader.UnitIdLoader()
        L.import_unit_ids(self.pdbs, recalculate=True)

    def test_import(self):
        ids = models.session.query(models.PdbUnitIdCorrespondence).\
                             filter_by(pdb='2Z4L').\
                             all()
        self.assertEqual(len(ids), 14256)

    def test_viral_import(self):
        ids = models.session.query(models.PdbUnitIdCorrespondence).\
                             filter_by(pdb='1A34').\
                             all()
        self.assertEqual(len(ids), 20557)

    def tearDown(self):
        models.session.query(models.PdbUnitIdCorrespondence).delete()


if __name__ == '__main__':
    unittest.main()