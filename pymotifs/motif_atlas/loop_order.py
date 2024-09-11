"""Load the loop ordering data.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.constants import NR_CACHE_NAME

from pymotifs.motif_atlas.utils import BaseLoader
from pymotifs.motif_atlas.info import Loader as InfoLoader
from pymotifs.motif_atlas.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, InfoLoader])
    @property
    def table(self):
        return mod.MlLoopOrder

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
