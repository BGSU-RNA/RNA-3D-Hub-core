"""Store the parents of nr classes.
"""

import datetime as dt

from pymotifs import core

from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.utils import BaseLoader
from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.id_mapping import Loader as IdMappingLoader


class Loader(BaseLoader):
    dependencies = set([IdMappingLoader, ChainLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days
    @property
    def table(self):
        return mod.NrClassParents

    def parents(self, release, grouping):
        data = []
        for group in grouping:
            for parent in group['parents']:
                data.append({
                    'nr_class_id': group['name']['class_id'],
                    'nr_release_id': release,
                    'nr_class_parent_id': parent['name']['class_id']
                })
        return data

    def no_parents(self, data):
        classes = [d['classes'] for d in data['parent_counts']]
        counts = [abs(d['unchanged']) + abs(d['updated']) for d in classes]
        return not sum(counts)

    def data(self, release, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("No grouping loaded")
        if data['parent'] == data['release']:
            raise core.Skip("First release has no parents")
        if self.no_parents(data):
            raise core.Skip("Parent counts shows no parents")
        return self.parents(release, data['groups'])
