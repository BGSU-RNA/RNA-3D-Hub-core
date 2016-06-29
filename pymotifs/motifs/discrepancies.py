"""This is a class to load all loop loop discrepancies that the motif atlas
produces. This relies upon there being cached data avaiable to store.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(core.SimpleLoader):
    dependencies = set([ReleaseLoader])
    table = mod.MlMutualDiscrepancy

    def query(self, session, *args, **kwargs):
        current, _ = ReleaseLoader(self.config, self.session).current_id()
        return session.query(mod.MlMutualDiscrepancy).\
            filter_by(release_id=current)

    def data(self):
        data = []
        for loop_type in ReleaseLoader.types:
            cached = self.cached(loop_type)
            data.extend(cached['discrepancies'])
        return data
