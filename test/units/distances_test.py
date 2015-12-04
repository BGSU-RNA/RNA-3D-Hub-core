from nose import SkipTest

from test import StageTest
from test import CifStageTest

import numpy as np

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

    def setUp(self):
        super(ComputingCentersTest, self).setUp()
        self.residues = list(self.structure.residues())

    def test_given_rna_gets_base_center(self):
        val = self.loader.center(self.residues[8])
        ans = np.array([-0.4305, -3.74425, 18.200875])
        np.testing.assert_array_almost_equal(ans, val)

    def test_given_aa_gets_backbone_center(self):
        raise SkipTest()

    def test_given_dna_gets_base_center(self):
        raise SkipTest()

    def test_given_anything_else_gets_overall_center(self):
        raise SkipTest()


class ComputingDistancesTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/124D.cif'

    def setUp(self):
        super(ComputingDistancesTest, self).setUp()
        self.residues = list(self.structure.residues())

    def test_can_compute_rna_rna_distance(self):
        val = self.loader.distance(self.residues[8], self.residues[9])
        ans = 4.266
        np.testing.assert_almost_equal(ans, val, decimal=1)

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
        val = list(self.loader.data(self.structure))
        self.assertEquals(3284, len(val))


class ProblematicStructureTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/1A34.cif'

    def test_can_load_distances(self):
        distances = list(self.loader.data(self.structure))
        self.assertTrue(distances)
