"""


"""


import unittest
import os


import PairwiseInteractionsLoader as loader
import models


class TestPairwiseInteractionsLoader(unittest.TestCase):

    def setUp(self):
        self.pdbs = ['1FG0','1HLX']
        self.loader = loader.PairwiseInteractionsLoader()

    def clean_up_database(self):
        """delete data from nr_tables"""
        session = models.session
        session.query(models.PairwiseInteractions).delete()
        session.commit()

    def test_import(self):
        self.loader.import_interactions(self.pdbs, recalculate=True)
        self.assertTrue( self.loader.success )


if __name__ == '__main__':
    unittest.main()
