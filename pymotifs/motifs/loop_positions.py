"""Load the loop to motif assignments.
"""

from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, InfoLoader])
    table = mod.MlLoopOrder

    def as_position(self, entry):
        return {
            'motif_id': entry['motif_id'],
            'loop_id': entry['positions']['loop_id'],
            'ml_release_id': entry['release_id'],
            'position': entry['positions']['position'],
            'unit_id': entry['positions']['unit_id'],
        }

    def data(self, pair, **kwargs):
        cached = self.cached(pair[0])
        return [self.as_position(entry) for entry in cached]
