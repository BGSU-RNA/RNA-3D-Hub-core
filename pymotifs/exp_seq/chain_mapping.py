from pymotifs import core

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.models import ExpSeqChainMapping as Mapping
from pymotifs.models import ChainInfo
from pymotifs.utils.structures import Structure

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.exp_seq.info import Loader as InfoLoader


class Loader(core.Loader):
    dependencies = set([ChainLoader, InfoLoader])

    def known(self, pdb):
        with self.session() as session:
            query = session.query(Mapping.chain_id).\
                join(ChainInfo, ChainInfo.chain_id == Mapping.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            return [result.chain_id for result in query]

    def remove(self, pdb, **kwargs):
        ids = self.known(pdb)
        if not ids:
            self.logger.info("Nothing to remove for %s", pdb)
            return None

        with self.session() as session:
            session.query(Mapping).\
                filter(Mapping.chain_id.in_(ids)).\
                delete(synchronize_session=False)

    def missing(self, pdb):
        helper = Structure(self.session.maker)
        possible = set(p[1] for p in helper.rna_chains(pdb, return_id=True))
        return possible - set(self.known(pdb))

    def has_data(self, pdb, **kwargs):
        return not self.missing(pdb)

    def exp_id(self, chain_id):
        with self.session() as session:
            query = session.query(Exp.exp_seq_id).\
                join(ChainInfo, ChainInfo.sequence == Exp.sequence).\
                filter(ChainInfo.chain_id == chain_id)
            return [result.id for result in query]

    def data(self, pdb, **kwargs):
        data = []
        missing = self.missing(pdb)
        for chain_id in missing:
            for exp_id in self.exp_id(chain_id):
                data.append(Mapping(exp_seq_id=exp_id, chain_id=chain_id))

        return data
