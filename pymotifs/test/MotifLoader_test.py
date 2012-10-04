"""

IL_20120905_0000: release 0.7, ilmarch_new
IL_20120906_1541: release 0.8 in MotifVersions_dev

"""


import unittest
import os
import logging


from MotifLoader import MotifLoader
import models


class TestMotifLoader(unittest.TestCase):

    ClassIsSetup = False

    def setUp(self):
        # If it was not setup yet, do it
        if not self.ClassIsSetup:
            print "Initializing testing environment"
            self.prepare()
            self.__class__.ClassIsSetup = True

    def prepare(self):
        unittest.TestCase.setUp(self)
        self.__class__.loader = MotifLoader(motif_type='IL')
        self.__class__.loader.start_logging()
        """override the config with the path to the test dataset"""
        script_path = os.path.dirname(os.path.abspath( __file__ ))
        self.__class__.loader.motifs_root = os.path.join(script_path,
                                                         'test_data',
                                                         'motifs')
        logging.info('Importing data from %s' % self.__class__.loader.motifs_root)
        """clean the database"""
        self.clean_up_database()
        """try importing the test dataset"""
        self.__class__.loader.import_data()

    def clean_up_database(self):
        """empty all motif-related tables"""
        session = models.session
        session.query(models.Release).delete()
        session.query(models.LoopOrder).delete()
        session.query(models.LoopPosition).delete()
        session.query(models.Motif).delete()
        session.query(models.Loop).delete()
        session.query(models.SetDiff).delete()
        session.query(models.Parents).delete()
        session.query(models.LoopDiscrepancy).delete()
        session.query(models.Release_diff).delete()
        session.commit()
        logging.info('Cleared old data from ml_tables')

    def destroy_environment(self):
        """
            remove all newly renamed mat files from all test release folders,
            delete text files with id correspondences
        """
        unittest.TestCase.tearDown(self)
        corr_file = 'correspondences.txt'
        # get the path to this script
        script_path   = os.path.dirname(os.path.abspath( __file__ ))
        # the test_data folder is in the same folder as this script
        releases_dir  = os.path.join(script_path, 'test_data', 'motifs')
        # loop over all releases
        for release_folder in os.listdir(releases_dir):
            release_folder = os.path.join(releases_dir, release_folder)
            if not os.path.isdir(release_folder):
                continue

            # remove new correspondences.txt files
            corr_file_path = os.path.join(release_folder, corr_file)
            if os.path.exists(corr_file_path):
                os.unlink(corr_file_path)
            else:
                logging.error('File %s not found' % corr_file_path)

            # remove all files in the mat folder
            mat_folder = os.path.join(release_folder, 'mat')
            for file in os.listdir(mat_folder):
                file_path = os.path.join(mat_folder, file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                except Exception, e:
                    logging.error(e)

    def test_import(self):
        """check that the import succeeded"""
        self.destroy_environment()
        self.assertTrue( self.loader.success )

    def test_release_count(self):
        """Release"""
        self.assertEqual(len(models.session.query(models.Release).all()), 2)

    def test_set_diff(self):
        """SetDiff"""
        self.assertEqual(len(models.session.\
                                    query(models.SetDiff).\
                                    filter(models.SetDiff.release_id=='0.2').\
                                    all()), 588)

    def test_releases_diff(self):
        """Release_diff"""
        self.assertEqual(len(models.session.\
                                    query(models.Release_diff).\
                                    all()), 1)

    def test_ml_loops(self):
        """Loop"""
        self.assertEqual(len(models.session.\
                                    query(models.Loop).\
                                    filter(models.Loop.release_id=='0.1').\
                                    all()), 1615)
        self.assertEqual(len(models.session.\
                                    query(models.Loop).\
                                    filter(models.Loop.release_id=='0.2').\
                                    all()), 1659)

    def test_ml_parents(self):
        """ml_history, Parents"""
        self.assertEqual(len(models.session.query(models.Parents).all()), 158)

    def test_ml_motifs(self):
        """Motif"""
        self.assertEqual(len(models.session.query(models.Motif).all()), 546)

    def test_mutual_disc_import(self):
        """ml_mutual_discrepancy, LoopDiscrepancy"""
        self.assertEqual(len(models.session.\
                                    query(models.LoopDiscrepancy).\
                                    all()), 192022)

    def test_loop_order_import(self):
        """ml_loop_order, LoopOrder"""
        self.assertEqual(len(models.session.\
                                    query(models.LoopOrder).\
                                    all()), 3274)


if __name__ == '__main__':
    unittest.main()