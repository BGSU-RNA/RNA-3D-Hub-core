from sqlalchemy.sql.expression import text

from pymotifs import core


ALL_PDBS_QUERY = """
select
    distinct P.pdb_id2
from correspondence_pdbs as P
join correspondence_info as I
on
    I.id = P.correspondence_id
where
    P.pdb_id1 = :pdb
    and I.good_alignment = 1
;
"""

ALL_CHAINS_QUERY = """
select
    P.*,
    C1.id as 'chain_id1',
    C2.id as 'chain_id2'
from correspondence_pdbs as P
join correspondence_info as I
on
    I.id = P.correspondence_id
join chain_info as C1
on
    C1.pdb_id = P.pdb_id1
    and C1.chain_name = P.chain_name1
join chain_info as C2
on
    C2.pdb_id = P.pdb_id2
    and C2.chain_name = P.chain_name2
where
    P.pdb_id1 = :pdb1
    and P.pdb_id2 = :pdb2
    and I.good_alignment = 1
    and C1.id != C2.id
;
"""


UNIT_ORDERING = """
select
    *
from correspondence_units
where
    pdb_id1 = :pdb1
    and pdb_id2 = :pdb2
    and chain1 = :chain1
    and chain2 = :chain2
    and correspondence_id = :corr_id
;
"""


class Helper(object):
    def __init__(self, maker):
        self.session = core.Session(maker)

    def pdbs(self, pdb):
        """Get all pdbs which have been aligned to the given pdb and whose
        alignment is good.
        """

        with self.session() as session:
            raw = text(ALL_PDBS_QUERY).bindparams(pdb=pdb)
            query = session.execute(raw).fetchall()
            data = [result.pdb_id2 for result in query]
            return data

    def chains(self, pdb1, pdb2):
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
            raw = text(ALL_CHAINS_QUERY).bindparams(pdb1=pdb1, pdb2=pdb2)
            query = session.execute(raw).fetchall()

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
            raw = text(UNIT_ORDERING).\
                bindparams(pdb1=chain1['pdb'], chain1=chain1['name'],
                           pdb2=chain2['pdb'], chain2=chain2['name'],
                           corr_id=corr_id)

            query = session.execute(raw)

            ordering = {}
            for result in query:
                ordering[result.unit_id1] = result.correspondence_index
                ordering[result.unit_id2] = result.correspondence_index
            return ordering
