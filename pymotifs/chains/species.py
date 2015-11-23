from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils.structures import Structure
from pymotifs.utils.structures import UnknownTaxonomyException

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.species_mapping import Loader as SpeciesLoader


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader, SpeciesLoader])
    table = mod.ChainSpecies

    def query(self, session, pdb):
        return session.query(mod.ChainSpecies).\
            join(mod.ChainInfo, mod.ChainInfo.chain_id == mod.ChainSpecies.chain_id).\
            filter(mod.ChainInfo.pdb_id == pdb)

    def validated(self, chain_id, assigned):
        given_name = None
        species = assigned
        with self.session() as session:
            given_name = session.query(mod.ChainInfo).get(chain_id).source

        if assigned in set([None, 512]) and given_name == 'Escherichia coli':
            self.logger.warning("Chain %s should probably have 562 as species",
                                chain_id)
            species = 562

        if species is None:
            return None
        return int(species)

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

            species = self.validated(chain_id, species)
            data.append({'chain_id': chain_id, 'species_id': species})
        return data
