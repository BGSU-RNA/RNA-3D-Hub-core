from pymotifs import core
from pymotifs.models import ChainSpecies
from pymotifs.models import ChainInfo
from pymotifs.utils.structures import Structure


class Loader(core.Loader):

    def has_data(self, pdb):
        with self.session() as session:
            query = session.query(ChainSpecies).\
                join(ChainInfo, ChainInfo.id == ChainSpecies.id).\
                filter(ChainInfo.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb):
        with self.session() as session:
            query = session.query(ChainInfo.id).filter(ChainInfo.pdb_id == pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            session.query(ChainSpecies).\
                filter(ChainSpecies.id.in_(ids)).\
                delete(synchronize_session=False)

    def data(self, pdb):
        helper = Structure(self.session.maker)
        data = []
        for chain_id, chain_name in helper.rna_chains(pdb, return_id=True):
            species = helper.source(pdb, chain_name, simplify=True)
            data.append(ChainSpecies(id=chain_id, species_id=species))
        return data
