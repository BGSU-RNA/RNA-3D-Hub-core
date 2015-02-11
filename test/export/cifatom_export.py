import os

from test import StageTest

from pymotifs.export.cifatom import Exporter


class ExporterFilenameTest(StageTest):
    loader_class = Exporter

    def test_can_generate_correct_pathname(self):
        val = self.loader.filename('03CX')
        ans = os.path.abspath("./FR3D/PDBFiles/03CX.cifatoms")
        self.assertEquals(ans, val)


class ExporterWritingTest(StageTest):
    loader_class = Exporter

    def has_cifatoms(self, pdb):
        return bool(self.loader.filename(pdb))

    def check_has_cifatoms(self, pdb):
        self.assertTrue(self.has_cifatoms(pdb))
        os.remove(self.loader.filename(pdb))

    def test_can_write_several_files(self):
        self.loader.process('1S72')
        self.check_has_cifatoms('1S72')
