"""This is a loader to load the chain to chain comparisons into the database.
This will look at good correspondences for a given structure and extract all
the aligned chains and then compute the geometric discrepancy between them and
then place them in the database.
"""

import itertools as it

from pymotifs import core
from pymotifs.utils import correspondence as corr

from pymotifs.models import ChainInfo
from pymotifs.models import ChainChainSimilarity as Similarity

from fr3d.geometry.discrepancy import discrepancy


class Loader(core.Loader):
    """A Loader to get all chain to chain similarity data.
    """
    allow_no_data = True

    def has_data(self, pdb, **kwargs):
        """Check if there are any chain_chain_similarity entries for this pdb.
        """

        with self.session() as session:
            query = session.query(Similarity).\
                join(ChainInfo, ChainInfo.id == Similarity.chain_id1).\
                filter(ChainInfo.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        """Remove all chain_chain_similarity entries for this pdb.
        """

        with self.session() as session:
            query = session.query(Similarity).\
                join(ChainInfo, ChainInfo.id == Similarity.chain_id1).\
                filter(ChainInfo.pdb_id == pdb)

            ids = [result.id for result in query]

        if not ids:
            self.logger.info("Nothing to remove for %s", pdb)
            return None

        with self.session() as session:
            session.query(Similarity).\
                filter(Similarity.id.in_(ids)).\
                delete(synchronize_session=False)

    def residues(self, cif, name, ordering):
        """Get the specified chain and extract only the residues in the
        ordering dictionary and return them in the specified order.

        :cif: The cif object to get the chain from.
        :chain: The name of the chain to access.
        :returns: A sorted list of the residues in the chain.
        """

        chain = cif.chain(1, name['name'])
        if not chain:
            raise core.InvalidState("Could not get chain %s for %s" %
                                    name, cif.pdb)

        residues = chain.residues()
        residues = it.ifilter(lambda r: r.unit_id() in ordering, residues)

        return sorted(residues, key=lambda r: ordering[r.unit_id()])

    def compare(self, corr_id, chain1, chain2):
        """Compare the chains. This will filter out all residues in the chain
        that do not have a base center computed.
        """

        residues1 = []
        residues2 = []
        for r1, r2 in zip(chain1['residues'], chain2['residues']):
            if 'base' in r1.centers and 'base' in r2.centers:
                residues1.append(r1)
                residues2.append(r2)

        if not residues1:
            self.logger.error("No residues with base centers for %s %s",
                              chain1['pdb'], chain1['name'])
            return None

        if not residues2:
            self.logger.error("No residues with base centers for %s %s",
                              chain2['pdb'], chain2['name'])
            return None

        self.logger.debug("Comparing %i residues in %s to %i residues in %s",
                          len(residues1), chain1['name'], len(residues2),
                          chain2['name'])

        disc = discrepancy(residues1, residues2)
        self.logger.debug("Got discrepancy %f", disc)

        return {
            'chain_id1': chain1['id'],
            'chain_id2': chain2['id'],
            'discrepancy': disc,
            'correspondence_id': corr_id
        }

    def data(self, pdb, **kwargs):
        """Compute all chain to chain similarity data. This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.
        """

        data = []

        util = corr.Helper(self.session.maker)
        cif1 = self.structure(pdb)
        cif1.infer_hydrogens()

        for pdb2 in util.pdbs(cif1.pdb):
            self.logger.debug("Comparing to %s", pdb2)

            cif2 = self.structure(pdb2)
            cif2.infer_hydrogens()

            corresponding = util.chains(cif1.pdb, cif2.pdb)
            for corr_id, chain1, chain2 in corresponding:
                self.logger.debug("Comparing chains %s %s",
                                  chain1['name'], chain2['name'])

                ordering = util.ordering(corr_id, chain1, chain2)
                if not ordering:
                    raise core.InvalidState("No ordering for %s, %s" %
                                            (chain1, chain2))

                chain1['residues'] = self.residues(cif1, chain1, ordering)
                chain2['residues'] = self.residues(cif2, chain2, ordering)

                compare = self.compare(corr_id, chain1, chain2)
                if compare:
                    data.append(Similarity(**compare))

        return data
