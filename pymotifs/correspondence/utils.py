import logging
import itertools as it
import collections as coll

import core

import utils as ut
from models import PdbCoordinates
from models import LoopsAll
from models import LoopPositions
from models import CorrespondenceInfo
from models import CorrespondenceNts


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

logger = logging.getLogger(__name__)


class StructureUtil(object):
    """Some useful methods
    """
    def __init__(self, maker):
        self.session = core.Session(maker)

    def longest_chain(self, pdb, model=1):
        with self.session() as session:
            query = session.query(PdbCoordinates).\
                filter_by(pdb=pdb, model=model).\
                filter(PdbCoordinates.unit.in_(['A', 'C', 'G', 'U'])).\
                order_by(PdbCoordinates.chain)

        grouped = it.groupby(query, lambda a: a.chain)
        max_pair = max(grouped, key=lambda (k, v): len(list(v)))
        return max_pair[0]

    def representative(self, pdb):
        with self.session() as session:
            result = session.execute(CURRENT_REP_QUERY, {'val': pdb})
            rep = set([rep[0] for rep in result.fetchall()])
        rep.discard(pdb)
        return rep

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

    def loops(self, pdb):
        """Get all loops in a structure.
        """
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

        if not mapping:
            logger.error("Could not generate mapping between %s to %s", ref,
                         pdb)

        return mapping
