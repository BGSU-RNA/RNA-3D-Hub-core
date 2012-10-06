"""

These tests make sure that Matlab clustering programs work. The results are
not imported into the database. That functionality is tested with a separate
test based on a test dataset.

"""


import unittest
import shutil
import os


import ClusterMotifs
import models


class TestClusterMotifs(unittest.TestCase):

    def setUp(self):
        """
            temporary copy over loop mat files from test_data into the
            standard location. Same for the file with a list of loops to cluster.
        """
        self.pdb_id = '1FG0'
        self.script_path = os.path.dirname(os.path.abspath( __file__ ))
        self.loader = ClusterMotifs.ClusterMotifs()
        self.loader.start_logging()
        self._copy_test_mat_files()

    def _copy_test_mat_files(self):
        # location of loop mat files
        src = os.path.join(self.script_path, 'test_data', 'clustering', self.pdb_id)
        # copy files to where Matlab expects them
        dest_folder = os.path.join(self.loader.config['locations']['loops_mat_files'], self.pdb_id)
        if not os.path.exists(dest_folder):
            os.makedirs(dest_folder)
        for filename in os.listdir(src):
            dest = os.path.join(dest_folder, filename)
            if not os.path.exists(dest):
                shutil.copyfile(os.path.join(src, filename), dest)

    def _list_loops(self, location, loop_type):
        """
            list all files in the location, return only loops of a certain type.
        """
        loop_ids = []
        files = os.listdir(location)
        for filename in files:
            if filename[:2] == loop_type:
                loop_ids.append(filename[:-4]) # discard .mat
        return loop_ids

    def _cluster_motifs(self, loop_type):
        """
        """
        self.loader.set_loop_type(loop_type)
        src = os.path.join(self.script_path, 'test_data', 'clustering', self.pdb_id)
        # manually set the list of loops to be clustered. Normally the list
        # comes from a database.
        self.loader.best_loops = self._list_loops(src, loop_type)
        self.loader.make_input_file_for_matlab()
        self.loader._make_release_directory()
        commands = self.loader.prepare_aAa_commands()
        self.loader.parallel_exec_commands( commands )
        self.loader.cluster_loops()

    def test_internal_loops(self):
        """
        """
        loop_type = 'IL'
        self._cluster_motifs(loop_type)
        self.assertTrue( self.loader.success )
        # compare lists of files

    def test_hairpin_loops(self):
        """
        """
        loop_type = 'HL'
        self._cluster_motifs(loop_type)
        self.assertTrue( self.loader.success )
        # compare lists of files


if __name__ == '__main__':
    unittest.main()
