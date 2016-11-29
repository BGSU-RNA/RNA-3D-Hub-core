from test import StageTest

from pymotifs import models as mod
from pymotifs.exp_seq.chain_mapping import Loader


class QueryingTest(StageTest):
    loader_class = Loader

    def chain_id(self, pdb, chain):
        with self.loader.session() as session:
            return session.query(mod.ChainInfo).\
                filter_by(pdb_id=pdb, chain_name=chain).\
                one().\
                chain_id

    def test_gets_correct_chains(self):
        assert self.loader.to_process(['1GID']) == sorted([
            self.chain_id('1GID', 'A'),
            self.chain_id('1GID', 'B'),
        ])


class DataTest(StageTest):
    loader_class = Loader

    def chain_id(self, pdb, chain):
        with self.loader.session() as session:
            return session.query(mod.ChainInfo).\
                filter_by(pdb_id=pdb, chain_name=chain).\
                one().\
                chain_id

    def test_can_get_correct_data(self):
        cid = self.chain_id('1GID', 'A')
        val = self.loader.data(cid)
        assert val.exp_seq_id == 40
        assert val.chain_id == cid
