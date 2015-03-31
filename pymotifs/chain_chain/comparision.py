"""This is a loader to load the chain to chain comparisons into the database.
This will look at good correspondences for a given structure and extract all
the aligned chains and then compute the geometric discrepancy between them.
"""

from sqlalchemy import aliased

from pymotifs import core

from pymotifs.models import ChainInfo
from pymotifs.models import CorrespondenceInfo as CorrInfo
from pymotifs.models import CorrespondencePdbs as CorrPdbs
from pymotifs.models import CorrespondenceUnits as CorrUnits
from pymotifs.models import ChainChainSimilarity as Similarity

from fr3d.geometry import discrepancy


class Loader(core.Loader):
    """A Loader to get all chain to chain similarity data.
    """
    allow_no_data = True

    def has_data(self, pdb, **kwargs):
        """Check if there are any chain_chain_similarity entries for this pdb.
        """

        with self.session() as session:
            query = session.query(Similarity).\
                join(ChainInfo, ChainInfo.id == Similarity.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        """Remove all chain_chain_similarity entries for this pdb.
        """

        with self.session() as session:
            query = session.query(Similarity).\
                join(ChainInfo, ChainInfo.id == Similarity.chain_id).\
                filter(ChainInfo.pdb_id == pdb)

            ids = [result.id for result in query]

        with self.session() as session:
            session.query(Similarity).\
                filter(Similarity.id.in_(ids)).\
                delete(synchronize_session=False)

    def corresponding_pdbs(self, pdb):
        """Get all pdbs which have been aligned to the given pdb and whose
        alignment is good.
        """

        with self.session() as session:
            query = session.query(CorrPdbs.pdb_id2).\
                join(CorrInfo, CorrInfo.id == CorrPdbs.correspondence_id).\
                filter_by(pdb_id1=pdb).\
                filter_by(good_alignment=True).\
                distinct()

            return [result.pdb_id2 for result in query]

    def corresponding_chains(self, pdb1, pdb2):
        """Get all chains which correspond between the two structures. This
        will return a list of 3 element tuples. The first will be the
        correspondence id, the second is a dict for chain1 and a second is a
        dict for chain2. Each chain dict will contain the name, the id, and the
        pdb.

        :params string pdb1: The first pdb.
        :params string pdb2: The second pdb.
        :returns: A list of tuples for the corresponding chains.
        """

        with self.session() as session:
            chain1 = aliased(ChainInfo)
            chain2 = aliased(ChainInfo)
            query = session.query(CorrPdbs,
                                  chain1.id.label('chain_id1'),
                                  chain2.id.label('chain_id2')).\
                filter_by(pdb_id1=pdb1, pdb_id2=pdb2)

            data = []
            for result in query:
                corr_id = result.correspondence_id
                chain1 = {'id': result.chain_id1, 'name': result.chain_name1}
                chain1['pdb'] = pdb1
                chain2 = {'id': result.chain_id2, 'name': result.chain_name2}
                chain2['pdb'] = pdb2
                data.append((corr_id, chain1, chain2))

        return data

    def ordering(self, corr_id, chain1, chain2):
        """Load the ordering of units in the given chain to chain
        correspondence. This will find the ordering of units in the
        correspondence from chain1 to chain2. The resulting dictionary will
        have entries for units in both chain1 and chain2. These entries may not
        start at 0 but are ordered to be increasing. Also they will not contain
        entries where either unit is not aligned.

        :param int corr_id: The correspondence id.
        :param dict chain1: The first chain.
        :param dict chain2: The second chain.
        :returns: An ordering dictionary.
        """

        with self.session() as session:
            query = session.query(CorrUnits).\
                filter_by(pdb_id1=chain1['pdb'], chain1=chain1['name']).\
                filter_by(pdb_id2=chain2['pdb'], chain2=chain2['name']).\
                filter_by(correspondence_id=corr_id)

            ordering = {}
            for result in query:
                ordering[result.unit_id1] = result.correspondence_index
                ordering[result.unit_id2] = result.correspondence_index
            return ordering

    def chain(self, cif, chain, ordering):
        """Get the specified chain and extract only the residues in the
        ordering dictionary and return them in the specified order.

        :cif: The cif object to get the chain from.
        :chain: The name of the chain to access.
        :returns: A sorted list of the residues in the chain.
        """

        residues = []
        for residue in cif.chain(chain).residues():
            if residue.unit_id in ordering:
                residues.append(residue)
        return sorted(residues, key=lambda r: ordering[r.unit_id])

    def data(self, pdb, **kwargs):
        """Compute all chain to chain similarity data. This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.
        """

        data = []

        cif1 = self.cif(pdb)
        for pdb2 in self.corresponding_pdbs(cif1.pdb):
            cif2 = self.cif(pdb2)

            corresponding = self.corresponding_chains(cif1.pdb, cif2.pdb)
            for corr_id, chain1, chain2 in corresponding:
                ordering = self.ordering(corr_id, chain1, chain2)
                chain1 = self.chain(cif1, chain1, ordering)
                chain2 = self.chain(cif2, chain2, ordering)

                data.append(Similarity(
                    chain_id1=chain1['id'],
                    chain_id2=chain2['id'],
                    discrepancy=discrepancy(chain1, chain2),
                    correspondence_id=corr_id
                ))

        return data
