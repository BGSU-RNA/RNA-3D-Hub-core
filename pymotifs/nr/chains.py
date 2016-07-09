import datetime as dt

from pymotifs import core

from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.utils import BaseLoader
from pymotifs.nr.id_mapping import Loader as IdMappingLoader


class Loader(BaseLoader):
    dependencies = set([IdMappingLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days
    table = mod.NrChains

    def chains(self, release, grouping):
        data = []
        for group in grouping:
            for chain in group['members']:
                data.append({
                    'ife_id': chain['id'],
                    'nr_class_id': group['name']['class_id'],
                    'nr_release_id': release,
                    'rank': chain['rank'],
                    'rep': chain == group['representative'],
                })
        return data

    def data(self, release, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("No grouping loaded")
        return self.chains(release, data['groups'])
