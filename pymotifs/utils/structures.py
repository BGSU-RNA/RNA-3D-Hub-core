"""Some queries that are useful in accessing stuff in the database.
"""

import itertools as it
import collections as coll

from sqlalchemy.orm import aliased

from pymotifs import core

from pymotifs import utils as ut
from pymotifs import models as mod

LONG_RANGE_CUTOFF = 4

POLYMER_UNITS_QUERY = '''
select *
from polymer_units as P
join pdb_unit_ordering as O
on
    O.nt_id = P.id
where
    P.pdb_id = :pdb
    and P.model = 1
    and P.chain = :chain
    and P.model = :model
order by P.polymer_id, O.index
'''

CURRENT_REP_QUERY = '''
select rep_id
from nr_current_representative
where
    pdb_id = :val
'''

CURRENT_MEMBERS_QUERY = '''
select N2.pdb_id
from nr_current_representative as N1
join nr_current_representative as N2
on
    N2.rep_id = N1.rep_id
where
    N1.pdb_id = :val
;
'''


class Base(object):
    def __init__(self, maker):
        self.session = core.Session(maker)


class Polymers(Base):
    """Methods for getting information about polymers.
    """

    def polymer_sequences(self, pdb, chain, model=1):
        results = self.__polymer_units__(pdb=pdb, chain=chain, model=model)
        sequence = []
        for _, nts in it.groupby(results, lambda r: r.polymer_id):
            sequence.append(''.join([nt.unit for nt in nts]))
        return sequence

    def unit_ids(self, pdb, chain, model=1):
        results = self.__polymer_units__(pdb=pdb, chain=chain, model=model)
        return [result.id for result in results]

    def __polymer_units__(self, **data):
        with self.session() as session:
            result = session.execute(POLYMER_UNITS_QUERY, data)

        Record = coll.namedtuple('Record', result.keys())
        return it.imap(lambda r: Record(*r), result.fetchall())


class NR(Base):
    def members(self, pdb):
        with self.session() as session:
            result = session.execute(CURRENT_MEMBERS_QUERY, {'val': pdb})
            rep = set([rep[0] for rep in result.fetchall()])
        rep.discard(pdb)
        return rep

    def representative(self, pdb):
        with self.session() as session:
            result = session.execute(CURRENT_REP_QUERY, {'val': pdb})
            rep = set([rep[0] for rep in result.fetchall()])
        rep.discard(pdb)
        return rep


class Structure(Base):
    def longest_chain(self, pdb, model=1):
        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter_by(pdb=pdb, model=model).\
                filter(mod.UnitInfo.unit.in_(['A', 'C', 'G', 'U'])).\
                order_by(mod.UnitInfo.chain)

        grouped = it.groupby(query, lambda a: a.chain)
        max_pair = max(grouped, key=lambda (k, v): len(list(v)))
        return max_pair[0]

    def interactions(self, pdb, chain):
        c1 = aliased(mod.UnitInfo)
        c2 = aliased(mod.UnitInfo)
        interactions = []
        with self.session() as session:
            query = session.query(mod.UnitPairsInteractions).\
                join(c1, c1.id == mod.UnitPairsInteractions.iPdbSig).\
                join(c2, c2.id == mod.UnitPairsInteractions.jPdbSig).\
                filter(mod.UnitPairsInteractions.pdb_id == pdb).\
                filter(c1.chain == c2.chain, c1.chain == chain)

            for result in query:
                data = ut.row2dict(result)
                data['id'] = int(data['id'])
                interactions.append(data)

        return interactions

    def stacks(self, pdb, chain, count=None):
        """Get the list or count of all stacking interactions.

        :pdb: The pdb file to query.
        :chain: The chain to get stacks for.
        :count
        """
        pass

    def loops(self, pdb):
        loops = []
        with self.session() as session:
            query = session.query(mod.LoopPositions).\
                join(mod.LoopsAll, mod.LoopPositions.loop_id == mod.LoopsAll.id).\
                filter(mod.LoopsAll.pdb == pdb).\
                order_by(mod.LoopsAll.id)

            grouped = it.groupby(it.imap(ut.row2dict, query),
                                 lambda a: a['loop_id'])
            for loop_id, positions in grouped:
                loops.append({
                    'id': loop_id,
                    'nts': [pos['nt_id'] for pos in positions]
                })

        return loops


class BasePairQueries(Base):
    """This is a class to deal with getting information about basepairs from
    the database. We store some useful but complex queries in here as methods
    for ease of use.
    """

    def representative(self, pdb, chain, count=False, range_cutoff=None,
                       near=False):
        """This gets all forward interactions within a chain. This means we get
        only the interactions which are either non-symmetric (tWH) or the
        symmetric ones where unit 1 id < unit 2 id.

        :pdb: The pdb to search.
        :chain: The chain to search in.
        :count: If we should return the count or the interactions.
        :returns: A count or the interactions themself.
        """
        inter = mod.UnitPairsInteractions

        with self.session() as session:
            u1, u2, query = self.__base__(session, pdb, chain, near=near)

            query = query.filter(u1.chain == u2.chain)

            if range_cutoff is not None:
                query = query.filter(inter.f_crossing >= range_cutoff)

            if count:
                return query.count()
            return [result for result in query]

    def cross_chain(self, pdb, chain, other_chain=None, count=False,
                    near=False):
        """This method gets all interactions between on chain and another
        chain. If no other chain is specified then we go to any other chain.
        This will get all interactions

        :pdb: The pdb to search.
        :chain: The first chain.
        :other_chain: The second chain(s) if any.
        :count: If we should return a count or all interactions
        :near: True if we should count nears or not
        :returns: The count or a list of bps.
        """

        with self.session() as session:
            u1, u2, query = self.__base__(session, pdb, chain, near=near)
            query = query.filter(u2.chain != u1.chain)

            if other_chain is not None:
                if isinstance(other_chain, list):
                    query = query.filter(u2.chain.in_(other_chain))
                else:
                    query = query.filter(u2.chain == other_chain)

            if count:
                return query.count()
            return [result for result in query]

    def between(self, pdb, chains, near=False, count=False):
        """This is a method to get the interactions between a list of chains.
        This is like a cross chain interaction but more general. Cross chain
        can be from A to B or C but not between B and C. This allows for
        interactions between A, B and C.

        :pdb: The pdb id.
        :chains: A list of chains.
        :near: A boolean to control showing nears.
        :count: A boolean to return a count or the list.
        :returns: A count or the list of interactions.
        """

        with self.session() as session:
            u1, u2, query = self.__base__(session, pdb, chains, near=near)
            query = query.\
                filter(u2.chain != u1.chain).\
                filter(u2.chain.in_(chains))

            if count:
                return query.count()
            return [result for result in query]

    def __base__(self, session, pdb, chain, symmetry=True, near=False):
        """A method to build the base queries for this class.

        :session: A database session.
        :pdb: The pdb id.
        :chain: The chain(s) to query.
        :symmetry: If we should care about symmetric base pairs.
        :near: If we should allow near interactions.
        :returns: A query.
        """

        u1 = aliased(mod.UnitInfo)
        u2 = aliased(mod.UnitInfo)
        bp = mod.BpFamilyInfo
        inter = mod.UnitPairsInteractions

        query = session.query(inter.unit1_id, inter.unit2_id, inter.f_lwbp).\
            join(u1, u1.id == inter.unit1_id).\
            join(u2, u2.id == inter.unit2_id).\
            join(bp, bp.id == inter.f_lwbp).\
            filter(bp.is_forward == True).\
            filter(inter.pdb_id == pdb)

        if symmetry:
            query = query.filter((bp.is_symmetric == False) |
                                 ((bp.is_symmetric == True) & (u1.id < u2.id)))

        if not near:
            query = query.filter(bp.is_near == False)

        if isinstance(chain, list):
            query = query.filter(u1.chain.in_(chain))
        else:
            query = query.filter(u1.chain == chain)

        return u1, u2, query


class Correspondence(Base):
    def reference(self, pdb):
        """Get all correlated reference structures.
        """
        with self.session() as session:
            query = session.query(mod.CorrespondenceInfo).filter_by(pdb2=pdb)
            return [ut.row2dict(result) for result in query]

    def mapping(self, ref, pdb):
        """Get the mapping from nucleotides in the reference to the nucleotides
        in the pdb.
        """
        mapping = {}
        with self.session() as session:
            query = session.query(mod.CorrespondenceNts).\
                join(mod.CorrespondenceInfo,
                     mod.CorrespondenceInfo.id ==
                     mod.CorrespondenceNts.correspondence_id).\
                filter(mod.CorrespondenceInfo.pdb1 == ref).\
                filter(mod.CorrespondenceInfo.pdb2 == pdb)
            for result in query:
                mapping[result.unit1_id] = result.unit2_id
                mapping[result.unit2_id] = result.unit1_id

        return mapping

    # def experimental_sequence_mapping(self, pdb, chain):
    #     with self.session() as session:
    #         query = session.query(ObsMap).\
    #             join(Obs, Obs.id == ObsMap.sequence_unit_id).\
    #             filter(ObsMap.pdb == pdb)
    #         return [ut.row2dict(result) for result in query]
