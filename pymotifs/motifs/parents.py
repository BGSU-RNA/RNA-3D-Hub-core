"""Load the motif parent data.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.motifs.release import Loader as ReleaseLoader
from pymotifs.constants import NR_CACHE_NAME


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, InfoLoader])
    table = mod.MlParents

    def parents(self, cached):
        data = []
        for motif in cached['motifs']:
            for parent in motif['parents']:
                data.append({
                    'ml_release_id': cached['release'],
                    'motif_id': motif['name']['full'],
                    'parent_ml_release_id': cached['parent'],
                    'parent_motif_id': parent['name']['full'],
                })
        return data

    def no_parents(self, data):
        classes = [d['motifs'] for d in data['parent_counts']]
        counts = [abs(d['unchanged']) + abs(d['updated']) for d in classes]
        return not sum(counts)

    def data(self, release, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("Missing cached data")

        if data['release'] == data['parent']:
            raise core.Skip("No parents for first release")
        if self.no_parents(data):
            raise core.Skip("Parent counts show no parents")

        return self.parents(data)
