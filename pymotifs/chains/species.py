from pymotifs import core
from pymotifs.models import ChainSpecies
from pymotifs.models import ChainInfo
from pymotifs.utils.structures import Structure
from pymotifs.utils.structures import UnknownTaxonomyException

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.species_mapping import Loader as SpeciesLoader


class Loader(core.Loader):
    dependencies = set([ChainLoader, SpeciesLoader])

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(ChainSpecies).\
                join(ChainInfo, ChainInfo.chain_id == ChainSpecies.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(ChainInfo.chain_id).\
                filter(ChainInfo.pdb_id == pdb)
            ids = [result.chain_id for result in query]

        with self.session() as session:
            session.query(ChainSpecies).\
                filter(ChainSpecies.chain_id.in_(ids)).\
                delete(synchronize_session=False)

    def data(self, pdb, **kwargs):
        helper = Structure(self.session.maker)
        data = []
        for chain_name, chain_id in helper.rna_chains(pdb, return_id=True):
            species = None
            try:
                species = helper.source(pdb, chain_name, simplify=True)
            except UnknownTaxonomyException as err:
                self.logger.warning(str(err))
            data.append(ChainSpecies(chain_id=chain_id, species_id=species))
        return data
