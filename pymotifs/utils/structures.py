"""Some queries that are useful in accessing stuff in the database.
"""

import itertools as it
import collections as coll

from sqlalchemy.orm import aliased

import core

import utils as ut
from models import PdbCoordinates
from models import LoopsAll
from models import LoopPositions
from models import CorrespondenceInfo
from models import CorrespondenceNts
from models import UnitPairsInteractions


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
            query = session.query(PdbCoordinates).\
                filter_by(pdb=pdb, model=model).\
                filter(PdbCoordinates.unit.in_(['A', 'C', 'G', 'U'])).\
                order_by(PdbCoordinates.chain)

        grouped = it.groupby(query, lambda a: a.chain)
        max_pair = max(grouped, key=lambda (k, v): len(list(v)))
        return max_pair[0]

    def interactions(self, pdb, chain):
        c1 = aliased(PdbCoordinates)
        c2 = aliased(PdbCoordinates)
        interactions = []
        with self.session() as session:
            query = session.query(UnitPairsInteractions).\
                join(c1, c1.id == UnitPairsInteractions.iPdbSig).\
                join(c2, c2.id == UnitPairsInteractions.jPdbSig).\
                filter(UnitPairsInteractions.pdb_id == pdb).\
                filter(c1.chain == c2.chain, c1.chain == chain)

            for result in query:
                data = ut.row2dict(result)
                data['id'] = int(data['id'])
                interactions.append(data)

        return interactions

    def loops(self, pdb):
        loops = []
        with self.session() as session:
            query = session.query(LoopPositions).\
                join(LoopsAll, LoopPositions.loop_id == LoopsAll.id).\
                filter(LoopsAll.pdb == pdb).\
                order_by(LoopsAll.id)

            grouped = it.groupby(it.imap(ut.row2dict, query),
                                 lambda a: a['loop_id'])
            for loop_id, positions in grouped:
                loops.append({
                    'id': loop_id,
                    'nts': [pos['nt_id'] for pos in positions]
                })

        return loops


class Correspondence(Base):
    def reference(self, pdb):
        """Get all correlated reference structures.
        """
        with self.session() as session:
            query = session.query(CorrespondenceInfo).filter_by(pdb2=pdb)
            return [ut.row2dict(result) for result in query]

    def mapping(self, ref, pdb):
        """Get the mapping from nucleotides in the reference to the nucleotides
        in the pdb.
        """
        mapping = {}
        with self.session() as session:
            query = session.query(CorrespondenceNts).\
                join(CorrespondenceInfo,
                     CorrespondenceInfo.id ==
                     CorrespondenceNts.correspondence_id).\
                filter(CorrespondenceInfo.pdb1 == ref).\
                filter(CorrespondenceInfo.pdb2 == pdb)
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
