import operator as op
import itertools as it
import collections as coll

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
                filter(info.good_alignment >= 1).\
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
                filter(corr.good_alignment >= 1).\
                filter(pdbs.chain_id_1 != pdbs.chain_id_2).\
                filter(ife1.is_integral == 1).\
                filter(ife2.is_integral == 1)

            data = []
            for result in query:
                chain1 = {
                    'pdb': pdb1,
                    'id': result.chain_id_1,
                    'name': result.chain_name_1
                }
                chain2 = {
                    'pdb': pdb2,
                    'id': result.chain_id_2,
                    'name': result.chain_name_2
                }
                data.append((result.correspondence_id, chain1, chain2))

        key = op.itemgetter('pdb', 'name')
        return sorted(data, key=lambda p: (p[0], key(p[1]), key(p[2])))

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

    def aligned_chains(self, ids, good=None, use_names=False):
        """
        Determine which chains a good alignment between them. This will
        produce a dictionary of dictionaries where the final values are a
        boolean indicating a good alignment or not.

        If good is given and not None then this will only load the alignments
        which have alignments marked good, or those bad. This means that the
        produced matrix will not contain all possible pairs of values. Specifically
        it will not contain pairs of chains that were never aligned when good
        is False. We will have entries at the top level for all given chains
        though, even if they have no matching alignments.

        :param list ids: A list of pdb ids to use.
        :param bool good: Bool to control if this should only load good
        alignments. True means only good alignments, False means only bad, and
        None means all alignments.
        :param bool use_names: A flag to indicate if we should use the names
        for chains or the id stored in the database.
        :returns: A dictionary of dictionaries indicating a good alignment or
        not.
        """

        with self.session() as session:
            exps = mod.ExpSeqPdb
            query = session.query(exps).\
                filter(exps.pdb_id.in_(ids))

            exp_mapping = coll.defaultdict(set)
            for result in query:
                name = result.chain_id
                if use_names:
                    name = '%s||%s' % (result.pdb_id, result.chain_name)
                exp_mapping[result.exp_seq_id].add(name)

        with self.session() as session:
            info = mod.CorrespondenceInfo
            query = session.query(info)

            if good is not None:
                if good:
                    query = query.filter(info.good_alignment >= 1)
                else:
                    query = query.filter(info.good_alignment == good)

            mapping = coll.defaultdict(dict)
            for result in query:
                is_good = bool(result.good_alignment)
                pairs = it.product(exp_mapping[result.exp_seq_id_1],
                                   exp_mapping[result.exp_seq_id_2])
                # pairs = it.ifilter(lambda (n1, n2): n1 != n2, pairs)
                pairs = it.ifilter(lambda pair: pair[0] != pair[1], pairs)
                for name1, name2 in pairs:
                    mapping[name1][name2] = is_good
                    mapping[name2][name1] = is_good

        if good is None:
            ids = it.chain.from_iterable(exp_mapping.values())
            pairs = it.product(ids, repeat=2)
            # pairs = it.ifilter(lambda (n1, n2): n1 != n2, pairs)
            pairs = it.ifilter(lambda pair: pair[0] != pair[1], pairs)
            pairs = list(pairs)
            for name1, name2 in pairs:
                if name2 not in mapping[name1]:
                    mapping[name1][name2] = False
                    mapping[name2][name1] = False
        else:
            for name in it.chain.from_iterable(exp_mapping.itervalues()):
                if name not in mapping:
                    mapping[name] = {}

        return dict(mapping)

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
