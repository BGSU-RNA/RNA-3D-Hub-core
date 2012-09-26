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
from nratlas import NRSqlAlchemyClasses


class TestNRClassLoader(unittest.TestCase):

    def setUp(self):
        """
            drop and recreate the NR database tables - should only affect the
            _test database as specified in the config file.
        """
        NRSqlAlchemyClasses.drop_all()
        NRSqlAlchemyClasses.create_all()
        # create object instance
        self.loader = NRClassLoader.Loader()
        # point to the test data directory
        script_path = os.path.dirname(os.path.abspath( __file__ ))
        self.loader.nrlists_root = os.path.join(script_path, 'test_data', 'nrlist')
        self.loader.start_logging()
        self.loader.import_data()

    def test_import(self):
        self.assertTrue( self.loader.success )

    def test_release_count(self):
        self.assertEqual(len(NRSqlAlchemyClasses.\
                                         session.\
                                         query(NRSqlAlchemyClasses.NR_release).\
                                         all()), 3)

    def test_set_diff(self):
        self.assertEqual(len(NRSqlAlchemyClasses.\
                                         session.\
                                         query(NRSqlAlchemyClasses.NR_setdiff).\
                                         all()), 46)

    def test_releases_diff(self):
        self.assertEqual(len(NRSqlAlchemyClasses.\
                                         session.\
                                         query(NRSqlAlchemyClasses.NR_release_diff).\
                                         all()), 24)

    def test_nr_pdbs(self):

        self.assertEqual(len(NRSqlAlchemyClasses.\
                                         session.\
                                         query(NRSqlAlchemyClasses.NR_pdb).\
                                         all()), 27484)

    def test_nr_parents(self):
        self.assertEqual(len(NRSqlAlchemyClasses.\
                                         session.\
                                         query(NRSqlAlchemyClasses.NR_parents).\
                                         all()), 23)

    def test_nr_classes(self):
        self.assertEqual(len(NRSqlAlchemyClasses.\
                                         session.\
                                         query(NRSqlAlchemyClasses.NR_class).\
                                         all()), 11067)


if __name__ == '__main__':
    unittest.main()
