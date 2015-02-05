from nose import SkipTest

from test import StageTest

from fr3d.cif.reader import Cif

from pymotifs.units.distances import Loader


class ComputingDistancesTest(StageTest):
    loader_class = Loader

    def test_can_compute_rna_rna_distance(self):
        raise SkipTest()

    def test_can_compute_rna_to_water_distance(self):
        raise SkipTest()

    def test_computes_rna_to_any_distance(self):
        raise SkipTest()

    def test_computes_protein_protein_distance(self):
        raise SkipTest()

    def test_computes_protein_rna_distance(self):
        raise SkipTest()

    def test_computes_water_water_distance(self):
        raise SkipTest()


class DistancesLoaderTest(StageTest):
    loader_class = Loader

    @classmethod
    def setUpClass(cls):
        with open('files/cif/1GID.cif', 'rb') as raw:
            cls.structure = Cif(raw).structure()

    # def setUp(self):
    #     super(DistancesLoaderTest, self).setUp()
    #     self.structure = self.__class__.structure
    #     self.distances = list(self.loader.data(self.structure))

    def test_can_load_distances(self):
        raise SkipTest()
        self.assertEquals(3284, len(self.distances))
