import os

from test import StageTest

from pymotifs.download import Downloader


class DownloadingTest(StageTest):
    loader_class = Downloader

    def test_can_get_cif_files(self):
        if os.path.exists("./FR3D/PDBFiles/1GID.cif"):
            os.remove("./FR3D/PDBFiles/1GID.cif")

        self.assertFalse(os.path.exists("./FR3D/PDBFiles/1GID.cif"))
        self.loader.process("1GID")
        self.assertTrue(os.path.exists("./FR3D/PDBFiles/1GID.cif"))

    def test_will_use_pdb_file_if_specified(self):
        ans = os.path.abspath('./FR3D/PDBFiles/1S72.pdb')
        self.assertEquals(ans, self.loader.filename('1S72', use_pdb=True))

    def test_it_will_use_pdb_url_if_specified(self):
        ans = 'http://www.rcsb.org/pdb/files/1S72.pdb.gz'
        self.assertEquals(ans, self.loader.url('1S72', use_pdb=True))
