"""Store the parents of nr classes.
"""

import datetime as dt

from pymotifs import core

from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.utils import BaseLoader
from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.class_rank import Loader as ClassRankLoader
from pymotifs.nr.id_mapping import Loader as IdMappingLoader


class Loader(BaseLoader):
    dependencies = set([IdMappingLoader, ClassRankLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days
    allow_no_data = True
    mark = False

    @property
    def table(self):
        return mod.NrClassParents

    def parents(self, release, grouping):
        data = []
        for group in grouping:
            #self.logger.info('Found group in nr.parents: %s', group)
            for parent in group['parents']:
                #self.logger.info('Found parent in nr.parents stage: %s', parent)
                if not parent:
                    self.logger.info('New group %s' % group['name']['full'])
                    data.append({
                        'nr_class_id': group['name']['class_id'],
                        'nr_release_id': release,
                        'nr_class_parent_id': parent['name']['class_id']
                    })
                else:
                    # self.logger.info('Exact match for %s, no need to generate new' % group['name']['full'])
                    pass
        return data

    def no_parents(self, data):
        classes = [d['classes'] for d in data['parent_counts']]
        counts = [abs(d['unchanged']) + abs(d['updated']) for d in classes]
        return not sum(counts)

    def data(self, release, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.Skip("Nothing to do here, maybe too few files, maybe an earlier stage failed")
        if data['parent'] == data['release']:
            raise core.Skip("First release has no parents")
        if self.no_parents(data):
            raise core.Skip("Parent counts shows no parents")
        return self.parents(release, data['groups'])
