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
