"""


"""


import unittest
import os
import logging


from DistancesAndCoordinatesLoader import DistancesAndCoordinatesLoader
from PdbInfoLoader import PdbInfoLoader
from LoopExtractor import LoopExtractor
from LoopQualityChecker import LoopQualityChecker
from PairwiseInteractionsLoader import PairwiseInteractionsLoader
from RedundantNucleotidesLoader import RedundantNucleotidesLoader
from BestChainsAndModelsLoader import BestChainsAndModelsLoader
from MotifAtlasBaseClass import MotifAtlasBaseClass
import models


class TestMotifAtlas(unittest.TestCase):

    def setUp(self):
        """runs the entire pipeline"""
        self.success = False

        self.clean_up_database()

        m = MotifAtlasBaseClass()
        os.chdir(m.config['locations']['fr3d_root'])
        m.start_logging()
        logging.info('Initializing update')

        """get new pdb files"""
        p = PdbInfoLoader()
        p.get_all_rna_pdbs()

        """override pdb files with a smaller set"""
        p.pdbs = ['1FG0','1HLX']

        """extract all loops and import into the database"""
        e = LoopExtractor()
        e.extract_and_import_loops(p.pdbs)

        """do loop QA, import into the database. Create a new loop release."""
        q = LoopQualityChecker()
        q.check_loop_quality(p.pdbs)

        """import pairwise interactions annotated by FR3D"""
        i = PairwiseInteractionsLoader()
        i.import_interactions(p.pdbs)

        """import coordinates and distances into the database"""
        d = DistancesAndCoordinatesLoader()
        d.import_distances(p.pdbs)
        d.import_coordinates(p.pdbs)

        """import info about redundant nucleotides"""
        r = RedundantNucleotidesLoader()
        r.import_redundant_nucleotides(p.pdbs)

        """import best chains and models"""
        b = BestChainsAndModelsLoader()
        b.import_best_chains_and_models(p.pdbs)

        self.success = True

    def test_import(self):
        self.assertTrue( self.success )

    def clean_up_database(self):
        """
        """
        session = models.session
        session.query(models.PdbAnalysisStatus).delete()
        session.query(models.LoopRelease).delete()
        session.query(models.LoopQA).delete()
        session.query(models.PairwiseInteractions).delete()
        session.query(models.Distances).delete()
        session.query(models.Coordinates).delete()
        session.query(models.AllLoops).delete()
        session.query(models.RedundantNucleotide).delete()
        session.query(models.PdbBestChainsAndModels).delete()
        session.commit()


    # test pdb_analysis_status dates.




if __name__ == '__main__':
    unittest.main()
