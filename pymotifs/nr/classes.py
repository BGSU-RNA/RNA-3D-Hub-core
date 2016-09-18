"""Store the nr classes (equivelance sets). This requires cached NR data to
process and then store.
"""

import datetime as dt

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.utils import BaseLoader
from pymotifs.nr.release import Loader as ReleaseLoader

from pymotifs.chains.loader import Loader as ChainLoader
from pymotifs.interactions.loader import Loader as InteractionLoader
from pymotifs.ife.loader import Loader as IfeLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader, ChainLoader, InteractionLoader,
                        IfeLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days

    @property
    def table(self):
        return mod.NrClasses

    def classes(self, data):
        result = []
        for group in data['groups']:
            result.append({
                'name': group['name']['full'],
                'nr_release_id': data['release'],
                'handle': group['name']['handle'],
                'version': group['name']['version'],
                'comment': group['comment'],
                'resolution': group['name']['cutoff']
            })
        return result

    def data(self, release_id, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("No cached data")
        return self.classes(data)
