from sqlalchemy.orm import aliased
from sqlalchemy.sql.expression import text

from pymotifs import core
from pymotifs import models as mod


ALL_PDBS_QUERY = """
select
    distinct P.pdb_id_2
from correspondence_pdbs as P
join correspondence_info as I
on
    I.correspondence_id = P.correspondence_id
where
    P.pdb_id_1 = :pdb
    and I.good_alignment = 1
;
"""

ALL_CHAINS_QUERY = """
select distinct
    P.correspondence_id,
    P.chain_name_1,
    P.chain_name_2,
    C1.id as 'chain_id1',
    C2.id as 'chain_id2'
from correspondence_pdbs as P
join correspondence_info as I
on
    I.correspondence_id = P.correspondence_id
join chain_info as C1
on
    C1.pdb_id = P.pdb_id_1
    and C1.chain_name = P.chain_name_1
join ife_chains as A1
on
    A1.chain_id = C1.id
join chain_info as C2
on
    C2.pdb_id = P.pdb_id_2
    and C2.chain_name = P.chain_name_2
join ife_chains as A2
on
    A2.chain_id = C2.id
where
    P.pdb_id_1 = :pdb1
    and P.pdb_id_2 = :pdb2
    and I.good_alignment = 1
    and C1.id != C2.id
    and A1.is_reference = 1
    and A2.is_reference = 1
;
"""


UNIT_ORDERING = """
select
    U.*
from correspondence_units as U
join chain_info as I1
on
    I1.pdb_id = U.pdb_id_1
    and I1.chain_name = U.chain_name_1
join chain_info as I2
on
    I2.pdb_id = U.pdb_id_2
    and I2.chain_name = U.chain_name_2
where
    I1.id = :chain1
    and I2.id = :chain2
    and correspondence_id = :corr_id
;
"""


UNIT_MAPPING = """
select
    *
from correspondence_units
where
    pdb_id_1 = :pdb1
    and pdb_id_2 = :pdb2
;
"""

ALIGNED_CHAINS = """
select distinct
    C1.chain_id as 'chain_id1',
    C2.chain_id as 'chain_id2',
    I.good_alignment as 'good_alignment'
from correspondence_pdbs as P
join correspondence_info as I
on
    I.corresponence_id = P.corresponence_id
join chain_info as C1
on
    C1.chain_id = P.chain_id_1
join chain_info as C2
on
    C2.chain_id = P.chain_id2
where
    P.pdb_id_1 in ({pdbs})
    and P.pdb_id_2 in ({pdbs})
    and C1.chain_id != C2.chain_id
;
"""


class Helper(core.Base):
    """A helper to load correspondence information. For example, getting what
    pdbs correspondence to other pdbs and such.
    """

    def pdbs(self, pdb):
        """Get all pdbs which have been aligned to the given pdb and whose
        alignment is good.
        """

        with self.session() as session:
            raw = text(ALL_PDBS_QUERY).bindparams(pdb=pdb)
            query = session.execute(raw).fetchall()
            data = [result.pdb_id_2 for result in query]
            return data

    def chains(self, pdb1, pdb2):
        """Get all chains which correspond between the two structures. This
        will return a list of 3 element tuples. The first will be the
        correspondence id, the second is a dict for chain1 and a second is a
        dict for chain2. Each chain dict will contain the name, the id, and the
        pdb. It will only load the reference chains from ife groups
        between the two structures.

        :params string pdb1: The first pdb.
        :params string pdb2: The second pdb.
        :returns: A list of tuples for the corresponding chains.
        """

        with self.session() as session:
            raw = text(ALL_CHAINS_QUERY).bindparams(pdb1=pdb1, pdb2=pdb2)
            query = session.execute(raw).fetchall()

            data = []
            for result in query:
                corr_id = result.id
                chain1 = {'id': result.chain_id_1, 'name': result.chain_name_1}
                chain1['pdb'] = pdb1
                chain2 = {'id': result.chain_id2, 'name': result.chain_name_2}
                chain2['pdb'] = pdb2
                data.append((corr_id, chain1, chain2))

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

    def pdb_mapping(self, pdb1, pdb2, reversed=False):
        """Get the mapping between two structures, ignoring the chain ids. This
        mapping is one way by default, that is pdb1 to pdb2 only, not the
        reverse. This can be changed with the reversed
        """

        with self.session() as session:
            raw = text(UNIT_MAPPING).\
                bindparams(pdb1=pdb1, pdb2=pdb2)
            results = session.execute(raw)

            mapping = {}
            corr_id = None
            for result in results:
                if corr_id is None:
                    corr_id = result.correspondence_id

                mapping[result.unit_id_1] = result.unit_id_2
                if reversed:
                    mapping[result.unit_id2] = mapping[result.unit_id_1]

            return corr_id, mapping

    def mapping(self, corr_id, chain1, chain2):
        """Get a mapping from nucleotides in chain1 to nucleotides in chain2.
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

        :ids: A list of pdb ids to use.
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
        with self.session() as session:
            raw = text(UNIT_ORDERING).\
                bindparams(chain1=chain_id1, chain2=chain_id2, corr_id=corr_id)

            query = session.execute(raw)
            for result in query:
                yield result
