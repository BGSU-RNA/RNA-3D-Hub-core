import pytest
from numpy.testing import assert_array_almost_equal

from test import StageTest

from pymotifs import core
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

    def test_knows_if_a_pair_has_been_done(self):
        assert self.loader.has_data(('1X8W|1|D', '1GRZ|1|B', 248L)) is True

    def test_knows_if_a_pair_is_not_done(self):
        assert self.loader.has_data(('1X8W|1|D', '1GRZ|1|B', -1)) is False

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

    def test_loads_an_ordering(self):
        val = self.loader.ordering(248, self.info1, self.info2)
        assert len(val) == 242 * 2
        assert val['1X8W|1|D|G|96'] == 0
        assert val['1GRZ|1|B|G|96'] == 0
        assert val['1X8W|1|D|G|414'] == 246
        assert val['1GRZ|1|B|G|414'] == 246

    def test_complains_if_unknown_corr_id(self):
        with pytest.raises(core.InvalidState):
            self.loader.ordering(-1, self.info1, self.info2)

    def test_it_loads_all_centers(self):
        ordering = self.loader.ordering(248, self.info1, self.info2)
        center1, center2 = self.loader.centers(248, self.info1, self.info2,
                                               ordering)
        assert len(center1) == 242
        assert len(center2) == 242
        assert_array_almost_equal(center1[0], [42.3341, 96.7131, 29.8361])
        assert_array_almost_equal(center2[0], [42.7088, 91.8778, 54.455])

    def test_loads_all_rotations(self):
        ordering = self.loader.ordering(248, self.info1, self.info2)
        rot1, rot2 = self.loader.rotations(248, self.info1, self.info2,
                                           ordering)
        assert len(rot1) == 242
        assert len(rot2) == 242
        assert_array_almost_equal(rot1[0], [[0.712491, -0.507095, 0.484986],
                                            [0.00272507, -0.689172, -0.724592],
                                            [0.701676, 0.517587, -0.489647]])
        assert_array_almost_equal(rot2[0], [[-0.0857475, -0.5474, 0.832466],
                                            [0.140197, -0.83386, -0.533876],
                                            [0.986404, 0.0709304, 0.148245]])

    def test_will_give_empty_for_unknown(self):
        ordering = self.loader.ordering(248, self.info1, self.info2)
        centers = self.loader.centers(None, self.info1, self.info2, ordering)
        rotations = self.loader.rotations(None, self.info1, self.info2,
                                          ordering)
        assert centers == ([], [])
        assert rotations == ([], [])


class ComputingDataTest(StageTest):
    loader_class = Loader

    def test_computes_both_discrepancies(self):
        val = self.loader.data(('1X8W|1|D', '1GRZ|1|B', 248L))
        self.assertTrue(len(val), 2)

    def test_it_increments_new_counts(self):
        assert self.loader.new_updates['1X8W|1|D'] == 0
        self.loader.data(('1X8W|1|D', '1GRZ|1|B', 248L))
        assert self.loader.new_updates['1X8W|1|D'] == 1
