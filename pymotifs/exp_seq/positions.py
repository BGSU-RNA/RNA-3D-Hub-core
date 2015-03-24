from pymotifs import core

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.models import ExpSeqPosition as Position
from pymotifs.models import ExpSeqChainMapping as Mapping
from pymotifs.models import ChainInfo


class Loader(core.Loader):

    def query(self, session, pdb):
        return session.query(Position.id).\
            join(Mapping, Position.exp_seq_id == Mapping.exp_seq_id).\
            join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
            filter(ChainInfo.pdb_id == pdb)

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = self.query(session, pdb)
            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = self.query(session, pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            query = session.query(Position.id).\
                filter(Position.id.in_(ids)).\
                delete(synchronize_session=False)

    def sequences(self, pdb):
        with self.session() as session:
            query = session.query(Exp).\
                join(Mapping, Exp.id == Mapping.exp_seq_id).\
                join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb)
            return [(result.id, result.sequence) for result in query]

    def positions(self, exp_id, sequence):
        positions = []
        for index, char in enumerate(sequence):
            positions.append({'exp_seq_id': exp_id,
                              'unit': char,
                              'index': index})
        return positions

    def data(self, pdb, **kwargs):
        data = []
        for (exp_id, sequence) in self.sequences(pdb):
            for position in self.positions(exp_id, sequence):
                data.append(Position(**position))
        return data
