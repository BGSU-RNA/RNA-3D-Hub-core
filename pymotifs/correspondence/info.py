from pymotifs import core
from pymotifs.models import ExpSeqInfo as ExpSeq
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.models import ExpSeqChainMapping as Mapping
from pymotifs.models import ChainInfo

import math


class MissingExpSeq(core.InvalidState):
    pass


class Loader(core.Loader):
    short_cutoff = 36

    def exp_seq(self, pdb):
        with self.session() as session:
            query = session.query(ExpSeq).\
                join(Mapping, Mapping.exp_seq_id == ExpSeq.id).\
                join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            sequences = []
            for result in query:
                sequences.append({
                    'id': result.id,
                    'length': result.length
                })

        return sequences

    def possible_short(self, exp_seq):
        with self.session() as session:
            query = session.query(Info).\
                filter(Info.length == exp_seq['length']).\
                filter(Info.id != exp_seq['id'])

            return [result.id for result in query]

    def possible_long(self, exp_seq):
        with self.session() as session:
            query = session.query(Info).\
                filter(Info.length <= 2 * exp_seq['length']).\
                filter(Info.length >= math.sqrt(exp_seq['length'])).\
                filter(Info.id != exp_seq['id'])

            return [result.id for result in query]

    def possible(self, exp_seq):
        if exp_seq['length'] < self.short_cutoff:
            return self.possible_short(exp_seq)
        return self.possible_long(exp_seq)

    def data(self, pdb):
        data = []
        for exp_seq in self.exp_seq(pdb):
            for other_seq in self.possible(exp_seq):
                data.append(Info(exp_seq_id1=exp_seq['id'],
                                 exp_seq_id2=other_seq))
        return data
