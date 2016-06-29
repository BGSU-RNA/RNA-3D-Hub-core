"""Load the loop ordering data.
"""

from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, InfoLoader])
    table = mod.MlLoopOrder

    def as_order(self, entry):
        return {
            'motif_id': entry['motif_id'],
            'loop_id': entry['ordering']['loop_id'],
            'ml_release_id': entry['release_id'],
            'orginal_order': entry['ordering']['orginal_order'],
            'similarity_order': entry['ordering']['similarity_order'],
        }

    def data(self, pair, **kwargs):
        cached = self.cached(pair[0])
        return [self.as_order(entry) for entry in cached]
