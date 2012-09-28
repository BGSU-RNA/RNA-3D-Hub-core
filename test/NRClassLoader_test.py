"""

the included nrlist test dataset corresponds to nrlists as follows:
20110528 - 0.19
20110604 - 0.20
20110611 - 0.12

the numbers in the assertEqual statements are given relative to the production
RNA 3D Hub database as of Sept 26, 2012.

"""


import unittest
import os


import NRClassLoader
import models


class TestNRClassLoader(unittest.TestCase):

    ClassIsSetup = False

    def setUp(self):
        # If it was not setup yet, do it
        if not self.ClassIsSetup:
            print "Initializing testing environment"
            self.prepare()
            self.__class__.ClassIsSetup = True

    def prepare(self):
        unittest.TestCase.setUp(self)
        self.clean_up_database()
        self.__class__.loader = NRClassLoader.Loader()
        # point to the test data directory
        script_path = os.path.dirname(os.path.abspath( __file__ ))
        self.__class__.loader.nrlists_root = os.path.join(script_path, 'test_data', 'nrlist')
        self.__class__.loader.start_logging()
        self.__class__.loader.import_data()

    def clean_up_database(self):
        """delete data from nr_tables"""
        session = models.session
        session.query(models.NR_release).delete(synchronize_session='fetch')
        session.query(models.NR_class).delete(synchronize_session='fetch')
        session.query(models.NR_pdb).delete(synchronize_session='fetch')
        session.query(models.NR_setdiff).delete(synchronize_session='fetch')
        session.query(models.NR_parents).delete(synchronize_session='fetch')
        session.query(models.NR_release_diff).delete(synchronize_session='fetch')
        session.query(models.NR_release_diff).delete(synchronize_session='fetch')
        session.commit()

    def test_import(self):
        self.assertTrue( self.loader.success )

    def test_release_count(self):
        self.assertEqual(len(models.session.query(models.NR_release).all()), 3)

    def test_set_diff(self):
        self.assertEqual(len(models.session.query(models.NR_setdiff).all()), 46)

    def test_releases_diff(self):
        self.assertEqual(len(models.session.query(models.NR_release_diff).all()), 24)

    def test_nr_pdbs(self):
        self.assertEqual(len(models.session.query(models.NR_pdb).all()), 27484)

    def test_nr_parents(self):
        self.assertEqual(len(models.session.query(models.NR_parents).all()), 23)

    def test_nr_classes(self):
        self.assertEqual(len(models.session.query(models.NR_class).all()), 11067)


if __name__ == '__main__':
    unittest.main()
