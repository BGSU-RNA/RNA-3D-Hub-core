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
    table = UnitInteractionSummary

    def query(self, session, pdb):
        return session.query(UnitInteractionSummary).\
            filter_by(pdb_id=pdb)

    def increment_bp(self, current, bp, crossing):
        return self.increment(current, 'bps', bp, crossing)

    def increment_stacks(self, current, stack, crossing):
        return self.increment(current, 'stacks', stack, crossing)

    def increment_bphs(self, current, unit1, unit2, bph, crossing):
        if unit1 != unit2 and bph != '0BPh':
            return self.increment(current, 'bphs', bph, crossing)
        return current

    def increment(self, current, family, name, crossing):
        if name and name[0] != 'n':
            current[name] += 1
            current['total'] += 1
            current[family] += 1
            lr_inc = int(crossing > LONG_RANGE)
            current['lr_' + name] += lr_inc
            current['lr_total'] += lr_inc
            current['lr_' + family] += lr_inc
        return current

    def data(self, pdb_id, **kwargs):

        with self.session() as session:
            query = session.query(UnitInfo.unit_id.label('unit_id_1'),
                                  UnitInfo.model,
                                  UnitInfo.chain,
                                  UnitInfo.pdb_id,
                                  UnitPairsInteractions.unit_id_2,
                                  UnitPairsInteractions.f_lwbp,
                                  UnitPairsInteractions.f_bphs,
                                  UnitPairsInteractions.f_stacks,
                                  UnitPairsInteractions.f_crossing,
                                  ).\
                outerjoin(UnitPairsInteractions,
                          UnitInfo.unit_id == UnitPairsInteractions.unit_id_1).\
                filter(UnitInfo.pdb_id == pdb_id).\
                filter(UnitInfo.unit_type_id == 'rna')

            data = coll.defaultdict(lambda: coll.defaultdict(int))
            for result in query:
                current = data[result.unit_id_1]
                current['unit_id'] = result.unit_id_1
                current['pdb_id'] = result.pdb_id
                current['model'] = result.model
                current['chain'] = result.chain
                crossing = result.f_crossing
                self.increment_bp(current, result.f_lwbp, crossing)
                self.increment_stacks(current, result.f_stacks, crossing)
                self.increment_bphs(current, result.unit_id_1,
                                    result.unit_id_2, result.f_bphs, crossing)
                data[current['unit_id']] = current

            return [(dict(v)) for v in data.values()]
