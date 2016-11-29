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

    def test_knows_if_file_exists(self):
        if not os.path.exists("./FR3D/PDBFiles/1GID.cif"):
            self.loader.process("1GID")
        assert self.loader.has_data('1GID') is True

    def test_knows_if_file_does_not_exist(self):
        assert self.loader.has_data('0GID') is False

    def test_can_remove_a_cif_file(self):
        if not os.path.exists("./FR3D/PDBFiles/1GID.cif"):
            self.loader.process("1GID")

        assert os.path.exists("./FR3D/PDBFiles/1GID.cif") is True
        self.loader.remove('1GID')
        assert os.path.exists("./FR3D/PDBFiles/1GID.cif") is False

    def test_will_get_cif_filename(self):
        ans = os.path.abspath('./FR3D/PDBFiles/1S72.cif')
        assert self.loader.filename('1S72') == ans

    def test_will_get_cif_url(self):
        assert self.loader.url('1S72') == 'http://www.rcsb.org/pdb/files/1S72.cif.gz'
