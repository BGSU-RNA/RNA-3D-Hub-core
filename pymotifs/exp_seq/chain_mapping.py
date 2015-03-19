from pymotifs import core

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.models import ExpSeqChainMapping as Mapping
from pymotifs.models import ChainInfo


class Loader(core.Loader):
    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(ChainInfo.id).filter_by(pdb_id=pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            session.query(Mapping).\
                filter(Mapping.chain_id.in_(ids)).\
                delete(synchronize_session=False)

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Mapping).\
                join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            return bool(query.count())

    def data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Exp.id, ChainInfo.id.label('chain_id')).\
                join(ChainInfo, ChainInfo.sequence == Exp.sequence).\
                filter(ChainInfo.pdb_id == pdb)

            data = []
            for result in query:
                data.append(Mapping(exp_seq_id=result.id,
                                    chain_id=result.chain_id))
            return data
