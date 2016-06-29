"""Load the motif info data.
"""

from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader])
    table = mod.MlMotifInfo

    def as_motif(self, entry):
        return self.table(
            motif_id=entry['motif_id'],
            ml_release_id=entry['release_id'],
            type=entry['loop_type'],
            handle=entry['name']['handle'],
            version=entry['name']['version'],
            comment=entry['name']['comment'],
        )

    def data(self, pair, **kwargs):
        cached = self.cached(pair[0])
        return [self.as_motif(entry) for entry in cached['motifs']]
