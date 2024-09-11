"""Load the loop to motif assignments.
"""

from pymotifs import models as mod

from pymotifs.motif_atlas.utils import BaseLoader
from pymotifs.motif_atlas.info import Loader as InfoLoader
from pymotifs.motif_atlas.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, InfoLoader])

    @property
    def table(self):
        return mod.MlLoopPositions

    def positions(self, cached):
        data = []
        for motif in cached['motifs']:
            for position in motif['positions']:
                data.append({
                    'motif_id': motif['motif_id'],
                    'loop_id': position['loop_id'],
                    'ml_release_id': cached['release'],
                    'position': position['position'],
                    'unit_id': position['unit_id'],
                })
        return data

    def data(self, pair, **kwargs):
        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("Missing cached data: %s", loop_type)
        return self.positions(cached)
