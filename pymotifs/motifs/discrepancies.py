"""This is a class to load all loop loop discrepancies that the motif atlas
produces. This relies upon there being cached data avaiable to store.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(core.SimpleLoader):
    dependencies = set([ReleaseLoader])
    table = mod.MlMutualDiscrepancy

    def to_process(self, *args, **kwargs):
        return ReleaseLoader.types

    def query(self, session, loop_type, **kwargs):
        data = self.cached(loop_type)
        return session.query(mod.MlMutualDiscrepancy).\
            join(mod.LoopInfo,
                 mod.LoopInfo.loop_id == mod.MlMutualDiscrepancy.loop_id_1).\
            filter(mod.LoopInfo.type == loop_type).\
            filter(mod.MlMutualDiscrepancy.ml_release_id == data['release'])

    def discrepancies(self, cached):
        data = []
        for disc in cached['discrepancies']:
            entry = dict(disc)
            entry['ml_release_id'] = cached['release']
            data.append(entry)
        return data

    def data(self, loop_type, **kwargs):
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No data cached")
        return self.discrepancies(cached)
