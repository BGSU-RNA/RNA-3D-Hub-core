"""Compute a summary of the number for each type of interaction for each unit.
This stage will compute the number of each type of basepair, base stack, and
base phosphate (with special cases for long range) that each unit in a
structure is the first entry of. This is useful in fr3d for searching as well
as for general querying and summarizing.
"""

import collections as coll

from pymotifs import core
from pymotifs.constants import LONG_RANGE

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsInteractions
from pymotifs.models import UnitInteractionSummary

from pymotifs.interactions.pairwise import Loader as InterLoader


class Loader(core.SimpleLoader):
    dependencies = set([InterLoader])

    def query(self, session, pdb):
        session.query(UnitInteractionSummary).\
            filter_by(pdb_id=pdb)

    def increment_bp(self, current, bp, crossing):
        self.increment(current, 'bps', bp, crossing)

    def increment_stacks(self, current, stack, crossing):
        self.increment(current, 'stacks', stack, crossing)

    def incremement_bphs(self, current, unit1, unit2, bph, crossing):
        if unit1 != unit2 and bph != '0BPh':
            self.increment(current, 'bphs', bph, crossing)

    def increment(self, current, family, name, crossing):
        if name and name[0] != 'n':
            current[name] += 1
            current['total'] += 1
            current[family] += 1
            if crossing > LONG_RANGE:
                current['lr_' + name] += 1
                current['lr_total'] += 1
                current['lr_' + family] += 1

    def data(self, pdb_id, **kwargs):

        with self.session() as session:
            query = session.query(UnitPairsInteractions,
                                  UnitInfo.model,
                                  UnitInfo.chain_name,
                                  ).\
                join(UnitInfo,
                     UnitInfo.unit_id == UnitPairsInteractions.unit_id_id).\
                filter_by(pdb_id=pdb_id)

            data = coll.defaultdict(lambda: coll.defaultdict(int))
            for result in query:
                current = result[result.unit_id_1]
                current['pdb_id'] = result.pdb_id
                current['model'] = result.model
                current['chain'] = result.chain_name
                crossing = result.f_crossing
                self.increment_bp(current, result.f_lwbp, crossing)
                self.increment_stacks(current, result.f_stacks, crossing)
                self.increment_bphs(current, result.f_bph, result.unit_id_1,
                                    result.unit_id_2, crossing)
                result[result.unit_id_1] = crossing

            return [UnitInteractionSummary(dict(v)) for v in data.values()]
