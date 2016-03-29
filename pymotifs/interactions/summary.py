import collections as coll

from pymotifs import core
from pymotifs import models as mod

from pymotifs.interactions.pairwise import Loader as PairwiseLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdb.info import Loader as PdbLoader


class Loader(core.SimpleLoader):
    dependencies = set([PdbLoader, UnitLoader, PairwiseLoader])

    table = mod.UnitInteractionSummary

    def query(self, session, pdb_id):
        return session.query(mod.UnitInteractionSummary).\
            filter_by(pdb_id=pdb)

    def interactions(self, pdb_id):
        with self.session() as session:
            query = session.query(mod.UnitPairsInteractions).\
                filter_by(pdb_id == pdb)

            interactions = []
            for result in query:
                is_lr = result.f_crossing >= 4
                annotations = [result.f_lwbp, result.f_stacks, result.f_bphs]
                annotations = [a for a in annotations if a]
                annotations = ['lr_' + a for a in a if is_lr]
                interactions.append({
                    'unit': result.unit_id_1,
                    'annotations': {
                        'bps': None,
                        'stacks': None,
                    }
                    'long_range': is_lr,
                })
        return interactions

    def data(self, pdb_id, **kwargs):
        summary = coll.defaultdict(lambda: defaultdict(int))
        for interaction in self.interactions(pdb_id):
            for annotation in interaction['annotations'].items():
                summary[unit1][annotation] += 1
                summary[unit1]['total'] += 1

        return [dict(summ) for summ in summary]
