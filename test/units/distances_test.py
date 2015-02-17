from nose import SkipTest

from test import StageTest
from test import CifStageTest

from pymotifs.units.distances import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_it_has_data(self):
        self.assertTrue(self.loader.has_data('124D'))

    def test_knows_has_no_data(self):
        self.assertFalse(self.loader.has_data('024D'))

    def test_can_delete_data(self):
        raise SkipTest()


class ComputingCentersTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/124D.cif'

    def test_given_rna_gets_base_center(self):
        raise SkipTest()

    def test_given_aa_gets_backbone_center(self):
        raise SkipTest()

    def test_given_dna_gets_base_center(self):
        raise SkipTest()

    def test_given_anything_else_gets_overall_center(self):
        raise SkipTest()


class ComputingDistancesTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/124D.cif'

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


class DistancesLoaderTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/1GID.cif'

    def test_can_load_distances(self):
        raise SkipTest()
        self.assertEquals(3284, len(self.distances))
