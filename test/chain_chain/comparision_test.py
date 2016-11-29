import pytest

from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal

from test import StageTest

from pymotifs import models as mod
from pymotifs.core import Skip
from pymotifs.utils import row2dict
from pymotifs.core import InvalidState
from pymotifs.chain_chain.comparision import Loader


class BaseTest(StageTest):
    loader_class = Loader

    def chain_id(self, pdb, chain):
        with self.loader.session() as session:
            return session.query(mod.ChainInfo).\
                filter_by(pdb_id=pdb, chain_name=chain).\
                one().\
                chain_id

    def pair(self, pdb1, chain1, pdb2, chain2):
        return (self.chain_id(pdb1, chain1), self.chain_id(pdb2, chain2))


class QueryTest(BaseTest):
    def test_it_knows_if_data_exists(self):
        val = self.pair('1ET4', 'A', '1ET4', 'B')
        assert self.loader.has_data(val) is True

    def test_it_knows_if_no_data_exists(self):
        assert self.loader.has_data((-1, 1)) is False

    def test_it_knows_if_reversed_data_exists(self):
        val = self.pair('1ET4', 'B', '1ET4', 'A')
        assert self.loader.has_data(val) is True


class ToProcessTest(BaseTest):
    @pytest.mark.skip()
    def test_it_creates_all_groups(self):
        pass

    @pytest.mark.skip()
    def test_it_filters_out_high_resolution(self):
        pass

    @pytest.mark.skip()
    def test_it_filters_out_no_resolution(self):
        pass


class LoadingInfoTest(BaseTest):
    def test_it_loads_correct_info(self):
        val = self.loader.info(4)
        assert val == {
            'pdb': '1CGM',
            'chain_id': self.chain_id('1CGM', 'I'),
            'chain_name': 'I',
            'model': 1,
            'sym_op': 'P_25',
            'ife_id': '1CGM|1|I',
            'name': '1CGM|1|I+P_25',
            'alt_id': None,
        }

    def test_it_complains_given_unknown_id(self):
        with pytest.raises(InvalidState):
            self.loader.info(-1)

    def test_complains_given_id_not_in_ife(self):
        with pytest.raises(InvalidState):
            self.loader.info(3)


class LoadingCorrespondneceTest(BaseTest):
    def test_can_load_a_correspondence(self):
        c1 = self.chain_id('4A3J', 'P')
        c2 = self.chain_id('3J9M', 'u')
        assert self.loader.corr_id(c1, c2) == 1

    def test_it_can_handle_reversed_ids(self):
        c1 = self.chain_id('4A3J', 'P')
        c2 = self.chain_id('3J9M', 'u')
        assert self.loader.corr_id(c2, c1) == 1

    def test_it_complains_given_invalid_pair(self):
        with pytest.raises(Skip):
            self.loader.corr_id(-1, 2)


class LoadingResiduesTest(BaseTest):
    def corr_id(self, info1, info2):
        return self.loader.corr_id(info1['chain_id'], info2['chain_id'])

    def setUp(self):
        super(LoadingResiduesTest, self).setUp()
        self.info1 = {
            'pdb': '1X8W',
            'chain_name': 'D',
            'chain_id': self.chain_id('1X8W', 'D'),
            'model': 1,
            'sym_op': '1_555',
            'alt_id': None,
        }
        self.info2 = {
            'pdb': '1GRZ',
            'chain_name': 'B',
            'chain_id': self.chain_id('1GRZ', 'B'),
            'model': 1,
            'sym_op': '1_555',
            'alt_id': None,
        }

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


class ComputingDataTest(BaseTest):
    def corr_id(self, chain1, chain2):
        return self.loader.corr_id(chain1, chain2)

    def test_computes_both_discrepancies(self):
        c1 = self.chain_id('1X8W', 'D')
        c2 = self.chain_id('1GRZ', 'B')
        corr_id = self.corr_id(c1, c2)
        val = [row2dict(d) for d in self.loader.data((c1, c2))]
        assert len(val) == 2
        # Remove discrepancy since it needs a different method
        assert_almost_equal(val[0].pop('discrepancy'), 0.227388, decimal=6)
        assert_almost_equal(val[1].pop('discrepancy'), 0.227388, decimal=6)
        del val[0]['chain_chain_similarity_id']
        del val[1]['chain_chain_similarity_id']
        assert val[0] == {
            'chain_id_1': c1,
            'chain_id_2': c2,
            'model_1': 1,
            'model_2': 1,
            'correspondence_id': corr_id,
            'num_nucleotides': 242
        }
        assert val[1] == {
            'chain_id_1': c2,
            'chain_id_2': c1,
            'model_1': 1,
            'model_2': 1,
            'correspondence_id': corr_id,
            'num_nucleotides': 242
        }
