import pytest

from test import StageTest
from test import CifStageTest

import numpy as np

from pymotifs.units.distances import Loader


class ToProcessTest(StageTest):
    loader_class = Loader

    def test_it_uses_all_given_pdbs_if_recalculate_is_true(self):
        pdbs = ['1S72', '1J5E', '2AW7']
        assert self.loader.to_process(pdbs, recalculate=True) == pdbs

    def test_it_uses_given_if_name_in_recalculate_list(self):
        pdbs = ['1S72', '1J5E', '2AW7']
        recalc = ['units.info', 'units.distances']
        assert self.loader.to_process(pdbs, recalculate=recalc) == pdbs

    def test_it_only_uses_new_pdbs_otherwise(self):
        pdbs = ['1S72', '1J5E', '1FJG', '4A3G', '1FEU', '5E3F', '2AW7']
        assert self.loader.to_process(pdbs, recalculate=False) == [
            '2AW7',
            '5E3F'
        ]

    def test_it_defaults_to_using_only_new(self):
        pdbs = ['1S72', '1J5E', '1FJG', '4A3G', '1FEU', '5E3F', '2AW7']
        assert self.loader.to_process(pdbs) == [
            '2AW7',
            '5E3F'
        ]


class QueryingTest(StageTest):
    loader_class = Loader

    def test_knows_it_has_data(self):
        self.assertTrue(self.loader.has_data('124D'))

    def test_knows_has_no_data(self):
        self.assertFalse(self.loader.has_data('024D'))

    @pytest.mark.skip()
    def test_can_delete_data(self):
        pass


class ComputingCentersTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/124D.cif'

    def test_given_rna_gets_base_center(self):
        print(list(self.structure.residues()))
        val = self.loader.center(self.structure.residue('124D|1|B|U|13'))
        ans = np.array([-5.94833333, 5.85976667, 8.00386667])
        np.testing.assert_array_almost_equal(ans, val)

    @pytest.mark.skip()
    def test_given_aa_gets_backbone_center(self):
        pass

    def test_given_dna_gets_center(self):
        val = self.loader.center(self.structure.residue('124D|1|A|DA|4'))
        ans = np.array([5.193938, 5.475406, 10.362031])
        np.testing.assert_array_almost_equal(ans, val)

    @pytest.mark.skip()
    def test_given_anything_else_gets_overall_center(self):
        pass


class ComputingDistancesTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/124D.cif'

    def dist(self, u1, u2):
        return self.loader.distance(self.structure.residue(u1),
                                    self.structure.residue(u2))

    def test_can_compute_rna_rna_distance(self):
        val = self.dist('124D|1|B|U|11', '124D|1|B|U|13')
        ans = 9.82
        np.testing.assert_almost_equal(ans, val, decimal=1)

    def test_computes_rna_to_any_distance(self):
        val = self.dist('124D|1|B|U|11', '124D|1|A|DA|4')
        ans = 15.07
        np.testing.assert_almost_equal(ans, val, decimal=1)

    @pytest.mark.skip()
    def test_computes_protein_protein_distance(self):
        pass

    @pytest.mark.skip()
    def test_computes_protein_rna_distance(self):
        pass


class DistancesLoaderTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/1GID.cif'

    def test_can_load_distances(self):
        val = list(self.loader.data(self.structure))
        self.assertEquals(2708, len(val))


class ProblematicStructureTest(CifStageTest):
    loader_class = Loader
    filename = 'test/files/cif/1A34.cif'

    def test_can_load_distances(self):
        distances = list(self.loader.data(self.structure))
        self.assertTrue(distances)
