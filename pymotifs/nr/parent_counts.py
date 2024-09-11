"""Import the parent counts.

This will compute the number of changes relative to the last release and put
then in the database.
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.utils import BaseLoader
from pymotifs.nr.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    dependencies = set([ReleaseLoader])

    mark = False

    @property
    def table(self):
        return mod.NrParentCounts

    def counts(self, release, parent, counts):
        data = []
        for count in counts:
            data.append({
                'resolution': count['cutoff'],
                'nr_release_id': release,
                'parent_nr_release_id': parent,
                'new_class_count': count['classes']['added'],
                'updated_class_count': count['classes']['updated'],
                'removed_class_count': count['classes']['removed'],
                'unchanged_class_count': count['classes']['unchanged'],
                'pdb_added_count': count['pdbs']['added'],
                'pdb_removed_count': count['pdbs']['removed'],
            })
        return data

    def data(self, release, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("No grouping loaded")
        return self.counts(release, data['parent'], data['parent_counts'])
