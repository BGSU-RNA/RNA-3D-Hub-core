"""Load the loop ordering data.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.constants import NR_CACHE_NAME

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

    def data(self, pair, **kwargs):
        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("Missing cached data: %s", loop_type)
        return self.ordering(cached)
