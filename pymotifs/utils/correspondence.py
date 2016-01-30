from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod


class Helper(core.Base):
    """A helper to load correspondence information. For example, getting what
    pdbs correspondence to other pdbs and such.
    """

    def pdbs(self, pdb):
        """Get all pdbs which have been aligned to the given pdb and whose
        alignment is good.

        :param str pdb: Pdb to get to get a list of all pdbs with corresponding
        chains.
        :returns: A list of pdb ids that have a corresponding chain to one
        chain in the given structure.
        """

        with self.session() as session:
            corr_pdb = mod.CorrespondencePdbs
            info = mod.CorrespondenceInfo
            query = session.query(corr_pdb.pdb_id_2).\
                join(info,
                     info.correspondence_id == corr_pdb.correspondence_id).\
                filter(corr_pdb.pdb_id_1 == pdb).\
                filter(info.good_alignment == 1).\
                order_by(corr_pdb.pdb_id_2).\
                distinct()
            return [result.pdb_id_2 for result in query]

    def chains(self, pdb1, pdb2):
        """Get all chains which correspond between the two structures. This
        will return a list of 3 element tuples. The first will be the
        correspondence id, the second is a dict for chain1 and a second is a
        dict for chain2. Each chain dict will contain the name, the id, and the
        pdb. It will only load the reference chains from ife groups
        between the two structures.

        :params str pdb1: The first pdb.
        :params str pdb2: The second pdb.
        :returns: A list of tuples for the corresponding chains.
        """

        with self.session() as session:
            pdbs = mod.CorrespondencePdbs
            corr = mod.CorrespondenceInfo
            ife1 = aliased(mod.IfeChains)
            ife2 = aliased(mod.IfeChains)
            query = session.query(pdbs.correspondence_id,
                                  pdbs.chain_name_1,
                                  pdbs.chain_name_2,
                                  pdbs.chain_id_1,
                                  pdbs.chain_id_2,
                                  ).\
                join(corr, corr.correspondence_id == pdbs.correspondence_id).\
                join(ife1, ife1.chain_id == pdbs.chain_id_1).\
                join(ife2, ife2.chain_id == pdbs.chain_id_2).\
                filter(pdbs.pdb_id_1 == pdb1).\
                filter(pdbs.pdb_id_2 == pdb2).\
                filter(corr.good_alignment == 1).\
                filter(pdbs.chain_id_1 != pdbs.chain_id_2).\
                filter(ife1.is_integral == 1).\
                filter(ife2.is_integral == 1)

            data = []
            for result in query:
                chain1 = {'id': result.chain_id_1, 'name': result.chain_name_1}
                chain1['pdb'] = pdb1
                chain2 = {'id': result.chain_id_2, 'name': result.chain_name_2}
                chain2['pdb'] = pdb2
                data.append((result.correspondence_id, chain1, chain2))

        return data

    def ordering(self, corr_id, chain_id1, chain_id2):
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

        ordering = {}
        for result in self.__ordering__(corr_id, chain_id1, chain_id2):
            ordering[result.unit_id_1] = result.correspondence_index
            ordering[result.unit_id_2] = result.correspondence_index
        return ordering

    def mapping(self, corr_id, chain1, chain2):
        """Get a mapping from nucleotides in chain1 to nucleotides in chain2.
        The mapping is in both directions chain1 to chain2 and chain2 to
        chain1.

        :param int corr_id: The correspondence id to use.
        :param dict chain1: A dictonary with an 'id' key for the chain id.
        :param dict chain2: A dictonary with an 'id' key for the chain id.
        :returns: A dictonary mapping between the units in the given chains.
        """

        mapping = {}
        for result in self.__ordering__(corr_id, chain1['id'], chain2['id']):
            mapping[result.unit_id_1] = result.unit_id_2
            mapping[result.unit_id_2] = result.unit_id_1
        return mapping

    def aligned_chains(self, ids, good=None):
        """Determine which chains a good alignment between them. This will
        produce a dictionary of dictionaries where the final values are a
        boolean indicating a good alignment or not.

        :param list ids: A list of pdb ids to use.
        :param bool good: Bool to control if this should only load good
        alignments.
        :returns: A dictionary of dictionaries indicating a good alignment or
        not.
        """

        mapping = {}
        with self.session() as session:
            c1 = aliased(mod.ChainInfo)
            c2 = aliased(mod.ChainInfo)
            m1 = aliased(mod.ExpSeqChainMapping)
            m2 = aliased(mod.ExpSeqChainMapping)
            info = mod.CorrespondenceInfo
            query = session.query(info.good_alignment,
                                  c1.chain_id.label('chain_id1'),
                                  c2.chain_id.label('chain_id2')).\
                join(m1, m1.exp_seq_id == info.exp_seq_id_1).\
                join(m2, m2.exp_seq_id == info.exp_seq_id_2).\
                join(c1, c1.chain_id == m1.chain_id).\
                join(c2, c2.chain_id == m2.chain_id).\
                filter(c1.pdb_id.in_(ids)).\
                filter(c2.pdb_id.in_(ids)).\
                filter(c1.chain_id != c2.chain_id)

            if good is not None:
                query = query.filter(info.good_alignment == int(good))

            for result in query:
                if result.chain_id1 not in mapping:
                    mapping[result.chain_id1] = {}
                if result.chain_id2 not in mapping:
                    mapping[result.chain_id2] = {}

                status = bool(result.good_alignment)
                mapping[result.chain_id1][result.chain_id2] = status
                mapping[result.chain_id2][result.chain_id1] = status

            return mapping

    def __ordering__(self, corr_id, chain_id1, chain_id2):
        """Compute the correspondence between the given chain ids and given
        correspondence id.

        :param int corr_id: The correpsondence id.
        :param int chain_id1: The first chain.
        :param int chain_id2: The second chain.
        :returns: A generator over the matching rows in correspondence_units.
        """

        with self.session() as session:
            units = mod.CorrespondenceUnits
            info1 = aliased(mod.ChainInfo)
            info2 = aliased(mod.ChainInfo)
            query = session.query(units).\
                join(info1, info1.pdb_id == units.pdb_id_1).\
                join(info2, info2.pdb_id == units.pdb_id_2).\
                filter(info1.chain_name == units.chain_name_1).\
                filter(info2.chain_name == units.chain_name_2).\
                filter(info1.chain_id == chain_id1).\
                filter(info2.chain_id == chain_id2).\
                filter(units.correspondence_id == corr_id)

            for result in query:
                yield result
