"""This is a class to load all loop loop discrepancies that the motif atlas
produces. This relies upon there being cached data avaiable to store.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motif_atlas.release import Loader as ReleaseLoader


class Loader(core.SimpleLoader):
    dependencies = set([ReleaseLoader])

    @property
    def table(self):
        return mod.MlMutualDiscrepancy

    def to_process(self, *args, **kwargs):
        if 'loop_type' in kwargs.get('manual', {}):
            loop_types = kwargs['manual']['loop_type'].split(',')
        else:
            loop_types = ReleaseLoader.types
        return loop_types

    def query(self, session, loop_type, **kwargs):
        data = self.cached(loop_type)
        if data:
            release = data['release']
        else:
            # early J8 releases, for example, had no data
            release = 'No release'

        return session.query(mod.MlMutualDiscrepancy).\
            join(mod.LoopInfo,
                 mod.LoopInfo.loop_id == mod.MlMutualDiscrepancy.loop_id_1).\
            filter(mod.LoopInfo.type == loop_type).\
            filter(mod.MlMutualDiscrepancy.ml_release_id == release)

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
            raise core.Skip("No %s data cached" % loop_type)
        discrepancy_data = self.discrepancies(cached)
        if discrepancy_data:
            return discrepancy_data
        else:
            raise core.Skip("No discrepancy data found, probably because the release is empty")
