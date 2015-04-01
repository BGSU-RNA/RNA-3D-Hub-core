from pymotifs import core

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.models import ExpSeqChainMapping as Mapping
from pymotifs.models import ChainInfo


class Loader(core.Loader):
    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Mapping.id).\
                join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb)
            ids = [result.id for result in query]

        if not ids:
            self.logger.info("Nothing to remove for %s", pdb)
            return None

        with self.session() as session:
            session.query(Mapping).\
                filter(Mapping.chain_id.in_(ids)).\
                delete(synchronize_session=False)

    def known(self, pdb):
        with self.session() as session:
            query = session.query(Mapping.chain_id).\
                join(ChainInfo, ChainInfo.id == Mapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            return [result.chain_id for result in query]

    def possible(self, pdb):
        with self.session() as session:
            query = session.query(ChainInfo.id).\
                join(Exp, Exp.sequence == ChainInfo.sequence).\
                filter(ChainInfo.pdb_id == pdb)
            return [result.id for result in query]

    def missing(self, pdb):
        return set(self.possible(pdb)) - set(self.known(pdb))

    def has_data(self, pdb, **kwargs):
        return not bool(self.missing(pdb))

    def data(self, pdb, **kwargs):
        data = []
        for chain_id in self.missing(pdb):
            with self.session() as session:
                query = session.query(Exp.id).\
                    join(ChainInfo, ChainInfo.sequence == Exp.sequence).\
                    filter(ChainInfo.id == chain_id)

                for result in query:
                    data.append(Mapping(exp_seq_id=result.id, chain_id=chain_id))

        return data
