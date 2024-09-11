"""Load the motif info data.

This will load the cached data to store all motifs into the DB.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motif_atlas.utils import BaseLoader
from pymotifs.motif_atlas.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader])
    @property
    def table(self):
        return mod.MlMotifsInfo

    def motifs(self, cached):
        data = []
        for entry in cached['motifs']:
            data.append(self.table(
                motif_id=entry['motif_id'],
                ml_release_id=cached['release'],
                type=cached['loop_type'],
                handle=entry['name']['handle'],
                version=entry['name']['version'],
                comment=entry['comment'],
            ))
        return data

    def data(self, pair, **kwargs):
        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No cached data")
        return self.motifs(cached)
