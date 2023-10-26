"""This stage will remove the cached nr data. This is run after all other
nr stages and only deletes the cached data.
"""

from pymotifs import core
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.classes import Loader as ClassLoader
# from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.parents import Loader as ParentLoader
from pymotifs.nr.parent_counts import Loader as CountLoader
# from pymotifs.nr.ordering import Loader as OrderingLoader
from pymotifs.nr.class_rank import Loader as ClassRankLoader


class Cleanup(core.MassLoader):
    allow_no_data = True
    dependencies = set([ReleaseLoader, ClassLoader, ParentLoader,
                        CountLoader, ClassRankLoader])

    def has_data(self, *args, **kwargs):
        grouping = self.cached(NR_CACHE_NAME)
        if not grouping:
            raise core.Skip("No precomputed grouping to cleanup")
        return False

    def remove(self, *args, **kwargs):
        self.logger.info("Nothing to remove")

    def data(self, *args, **kwargs):
        self.logger.info('Cleaning up nr data')
        self.evict(NR_CACHE_NAME)
        return None
