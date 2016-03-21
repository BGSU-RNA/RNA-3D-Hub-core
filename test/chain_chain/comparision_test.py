import pytest

from test import StageTest

from pymotifs.chain_chain.comparision import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def test_can_find_all_chains_to_process(self):
        val = self.loader.to_process(['1X8W', '1GRZ'])
        ans = [('1X8W|1|A', '1GRZ|1|A', 248L),
               ('1X8W|1|A', '1GRZ|1|B', 248L),
               ('1X8W|1|B', '1GRZ|1|A', 248L),
               ('1X8W|1|B', '1GRZ|1|B', 248L),
               ('1X8W|1|C', '1GRZ|1|A', 248L),
               ('1X8W|1|C', '1GRZ|1|B', 248L),
               ('1X8W|1|D', '1GRZ|1|A', 248L),
               ('1X8W|1|D', '1GRZ|1|B', 248L)]
        self.assertEquals(val, ans)

    @pytest.mark.skip()
    def test_knows_if_a_pair_has_been_done(self):
        self.assertFalse(self.loader.has_data(('1X8W|1|D', '1GRZ|1|B', 248L)))

    def test_knows_if_a_pair_is_not_done(self):
        self.assertFalse(self.loader.has_data(('1X8W|1|D', '1GRZ|1|B', 0)))


class LoadingResiduesTest(StageTest):
    loader_class = Loader

    @pytest.mark.skip()
    def test_loads_residues_in_order(self):
        pass

    @pytest.mark.skip()
    def test_complains_if_missing_residue(self):
        pass


class ComputingDataTest(StageTest):
    loader_class = Loader

    def test_computes_both_discrepancies(self):
        val = self.loader.data(('1X8W|1|D', '1GRZ|1|B', 248L))
        self.assertTrue(len(val), 2)
