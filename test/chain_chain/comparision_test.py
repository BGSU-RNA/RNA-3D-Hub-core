import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal

from test import StageTest

from pymotifs import models as mod
from pymotifs.chain_chain.comparision import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def corr_id(self, info1, info2):
        with self.loader.session() as session:
            pdbs = mod.CorrespondencePdbs
            return session.query(pdbs.correspondence_id).\
                filter(pdbs.pdb_id_1 == info1['pdb']).\
                filter(pdbs.pdb_id_2 == info2['pdb']).\
                filter(pdbs.chain_name_1 == info1['chain_name']).\
                filter(pdbs.chain_name_2 == info2['chain_name']).\
                one().\
                correspondence_id

    def test_can_find_all_chains_to_process(self):
        pairs = self.loader.to_process(['1X8W', '1GRZ'])
        val = [(p[0], p[1]) for p in pairs]
        ans = [('1GRZ|1|A', '1GRZ|1|B'),
               ('1X8W|1|A', '1GRZ|1|A'),
               ('1X8W|1|A', '1GRZ|1|B'),
               ('1X8W|1|A', '1X8W|1|B'),
               ('1X8W|1|A', '1X8W|1|C'),
               ('1X8W|1|A', '1X8W|1|D'),
               ('1X8W|1|B', '1GRZ|1|A'),
               ('1X8W|1|B', '1GRZ|1|B'),
               ('1X8W|1|B', '1X8W|1|C'),
               ('1X8W|1|B', '1X8W|1|D'),
               ('1X8W|1|C', '1GRZ|1|A'),
               ('1X8W|1|C', '1GRZ|1|B'),
               ('1X8W|1|C', '1X8W|1|D'),
               ('1X8W|1|D', '1GRZ|1|A'),
               ('1X8W|1|D', '1GRZ|1|B')]
        self.assertEquals(val, ans)

    def test_knows_if_a_pair_has_been_done(self):
        corr_id = self.corr_id({'pdb': '1X8W', 'chain_name': 'D'},
                               {'pdb': '1GRZ', 'chain_name': 'B'})
        assert self.loader.has_data(('1X8W|1|D', '1GRZ|1|B', corr_id)) is True

    def test_knows_if_a_pair_is_not_done(self):
        assert self.loader.has_data(('4V7R|1|A1', '4V88|1|A2', -1)) is False

    def test_it_will_have_data_if_done_many_new(self):
        assert self.loader.has_data(('1X8W|1|D', '1GRZ|1|B', 0)) is False
        self.loader.new_updates['1X8W|1|D'] = 11
        assert self.loader.has_data(('1X8W|1|D', '1GRZ|1|B', 0)) is True


class LoadingInfoTest(StageTest):
    loader_class = Loader

    def test_it_loads_all_info(self):
        val = self.loader.info('1X8W|1|D')
        ans = {
            'pdb': '1X8W',
            'chain_name': 'D',
            'chain_id': 251,
            'model': 1,
            'sym_op': '1_555',
        }
        assert val == ans


class LoadingResiduesTest(StageTest):
    loader_class = Loader

    def setUp(self):
        super(LoadingResiduesTest, self).setUp()
        self.info1 = self.loader.info('1X8W|1|D')
        self.info2 = self.loader.info('1GRZ|1|B')

    def corr_id(self, info1, info2):
        with self.loader.session() as session:
            pdbs = mod.CorrespondencePdbs
            return session.query(pdbs.correspondence_id).\
                filter(pdbs.pdb_id_1 == info1['pdb']).\
                filter(pdbs.pdb_id_2 == info2['pdb']).\
                filter(pdbs.chain_name_1 == info1['chain_name']).\
                filter(pdbs.chain_name_2 == info2['chain_name']).\
                one().\
                correspondence_id

    def test_it_loads_all_data(self):
        corr_id = self.corr_id(self.info1, self.info2)
        c1, c2, r1, r2 = self.loader.matrices(corr_id, self.info1, self.info2)
        assert len(c1) == 242
        assert len(c2) == 242
        assert len(r1) == 242
        assert len(r2) == 242
        assert_array_almost_equal(c1[0], [42.3341, 96.7131, 29.8361])
        assert_array_almost_equal(c2[0], [42.7088, 91.8778, 54.455])
        assert_array_almost_equal(r1[0], [[0.712491, -0.507095, 0.484986],
                                          [0.00272507, -0.689172, -0.724592],
                                          [0.701676, 0.517587, -0.489647]])
        assert_array_almost_equal(r2[0], [[-0.0857475, -0.5474, 0.832466],
                                          [0.140197, -0.83386, -0.533876],
                                          [0.986404, 0.0709304, 0.148245]])

    def test_will_give_empty_for_unknown(self):
        c1, c2, r1, r2 = self.loader.matrices(None, self.info1, self.info2)
        assert_array_almost_equal(c1, np.array([]))
        assert_array_almost_equal(c2, np.array([]))
        assert_array_almost_equal(r1, np.array([]))
        assert_array_almost_equal(r2, np.array([]))


class ComputingDataTest(StageTest):
    loader_class = Loader

    def corr_id(self, info1, info2):
        with self.loader.session() as session:
            pdbs = mod.CorrespondencePdbs
            return session.query(pdbs.correspondence_id).\
                filter(pdbs.pdb_id_1 == info1['pdb']).\
                filter(pdbs.pdb_id_2 == info2['pdb']).\
                filter(pdbs.chain_name_1 == info1['chain_name']).\
                filter(pdbs.chain_name_2 == info2['chain_name']).\
                one().\
                correspondence_id

    def test_computes_both_discrepancies(self):
        corr_id = self.corr_id({'pdb': '1X8W', 'chain_name': 'D'},
                               {'pdb': '1GRZ', 'chain_name': 'B'})
        val = self.loader.data(('1X8W|1|D', '1GRZ|1|B', corr_id))
        self.assertTrue(len(val), 2)

    def test_it_increments_new_counts(self):
        corr_id = self.corr_id({'pdb': '1X8W', 'chain_name': 'D'},
                               {'pdb': '1GRZ', 'chain_name': 'B'})
        assert self.loader.new_updates['1X8W|1|D'] == 0
        self.loader.data(('1X8W|1|D', '1GRZ|1|B', corr_id))
        assert self.loader.new_updates['1X8W|1|D'] == 1

    def test_it_can_handle_nan_discrepancy(self):
        corr_id = self.corr_id({'pdb': '1FCW', 'chain_name': 'A'},
                               {'pdb': '1FCW', 'chain_name': 'B'})
        val = self.loader.data(('1FCW|1|A', '1FCW|1|B', corr_id))
        assert_almost_equal(7.266e-05, val[0].discrepancy)
