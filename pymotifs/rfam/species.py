"""Assign species to each chain.

This will figure out what chain belongs to which species. This will map from
the assign taxon id to species id. It will also attempt to fix some mistakes
that authors make, such as labeling e. coli as 512 instead of 562. This can be
extended to deal with differences between taxon ids and organism names.
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils.structures import Structure
from pymotifs.utils.structures import UnknownTaxonomyException

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.species_mapping import Loader as SpeciesLoader


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader, SpeciesLoader])
    @property
    def table(self):
        return mod.ChainSpecies

    def query(self, session, pdb):
        return session.query(mod.ChainSpecies).\
            join(mod.ChainInfo, mod.ChainInfo.chain_id == mod.ChainSpecies.chain_id).\
            filter(mod.ChainInfo.pdb_id == pdb)

    def validated(self, chain_id, assigned):
        """Attempt to fix the assigned taxon id for a given chain id. This will
        look at data, like the name to produce a taxon id. For example things
        named 'Escherichia coli' and given label 512 (blue tongue virus) are
        labeled 562 (E coli) instead.

        :param int chain_id: Id of the chain to correct.
        :param int assigned: The currently assigned taxon id.
        :returns: A new taxon id.
        """

        given_name = None
        species = assigned
        with self.session() as session:
            given_name = session.query(mod.ChainInfo).get(chain_id).source

        if assigned in set([None, 512]) and given_name == 'Escherichia coli':
            self.logger.warning("Chain %s should probably have 562 as species" % chain_id)
            species = 562

        if species is None:
            return None
        return int(species)

    def data(self, pdb, **kwargs):
        """Compute and assignment of the chain ids in this pdb.

        :param str pdb: The PDB id.
        :returns: A dict with chain_id and species_id.
        """

        helper = Structure(self.session.maker)
        data = []
        rna_chains = helper.rna_chains(pdb, return_id=True)

        if not rna_chains:
            raise core.InvalidState("Structure %s contains no RNA" % pdb)

        for chain_name, chain_id in rna_chains:
            species = None
            try:
                species = helper.source(pdb, chain_name, simplify=True)
            except UnknownTaxonomyException as err:
                self.logger.warning(str(err))

            species = self.validated(chain_id, species)
            data.append({'chain_id': chain_id, 'species_id': species})
        return data
