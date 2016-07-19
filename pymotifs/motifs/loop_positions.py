"""Load the loop to motif assignments.
"""

from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, InfoLoader])
    table = mod.MlLoopPositions

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

    def data(self, release, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("Missing cached data")

        return self.positions(data)
