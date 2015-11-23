from pymotifs import core
from pymotifs.models import ChainSpecies
from pymotifs.models import ChainInfo
from pymotifs.utils.structures import Structure
from pymotifs.utils.structures import UnknownTaxonomyException

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.species_mapping import Loader as SpeciesLoader


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader, SpeciesLoader])

    def query(self, session, pdb):
        return session.query(ChainSpecies).\
            join(ChainInfo, ChainInfo.chain_id == ChainSpecies.chain_id).\
            filter(ChainInfo.pdb_id == pdb)

    def data(self, pdb, **kwargs):
        helper = Structure(self.session.maker)
        data = []
        rna_chains = helper.rna_chains(pdb, return_id=True)

        if not rna_chains:
            raise core.InvalidState("Structure %s contains no rna" % pdb)

        for chain_name, chain_id in rna_chains:
            species = None
            try:
                species = helper.source(pdb, chain_name, simplify=True)
            except UnknownTaxonomyException as err:
                self.logger.warning(str(err))

            species = species or 32630
            data.append(ChainSpecies(chain_id=chain_id, species_id=species))
        return data
