from pymotifs import core
from pymotifs.models import ChainSpecies
from pymotifs.models import ChainInfo
from pymotifs.utils.structures import Structure
from pymotifs.utils.structures import SYNTHEIC


class Loader(core.Loader):

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(ChainSpecies).\
                join(ChainInfo, ChainInfo.id == ChainSpecies.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(ChainInfo.id).filter(ChainInfo.pdb_id == pdb)
            ids = [result.chain_id for result in query]

        with self.session() as session:
            session.query(ChainSpecies).\
                filter(ChainSpecies.chain_id.in_(ids)).\
                delete(synchronize_session=False)

    def data(self, pdb, **kwargs):
        helper = Structure(self.session.maker)
        data = []
        for chain_name, chain_id in helper.rna_chains(pdb, return_id=True):
            species = helper.source(pdb, chain_name, simplify=True)
            if species is None:
                species = SYNTHEIC[0]
            data.append(ChainSpecies(chain_id=chain_id, species_id=species))
        return data
