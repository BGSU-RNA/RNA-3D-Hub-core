"""Load the loop ordering data.
"""

from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, InfoLoader])
    table = mod.MlLoopOrder

    def ordering(self, cached):
        data = []
        for entry in cached['motifs']:
            for loop in entry['ordering']:
                data.append({
                    'motif_id': entry['motif_id'],
                    'loop_id': loop['loop_id'],
                    'ml_release_id': cached['release'],
                    'original_order': loop['original_order'],
                    'similarity_order': loop['similarity_order'],
                })
        return data

    def data(self, release, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("Missing cached data")
        return self.ordering(cached)
