"""This stage will remove the cached nr data. This is run after all other
nr stages and only deletes the cached data.
"""

from pymotifs import core

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.classes import Loader as ClassLoader
from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.parents import Loader as ParentLoader
from pymotifs.nr.parent_counts import Loader as CountLoader


class Cleanup(core.Loader):
    allow_no_data = True
    dependencies = set([ReleaseLoader, ChainLoader, ClassLoader, ParentLoader,
                        CountLoader])

    def has_data(self, *args, **kwargs):
        grouping = self.cached('nr')
        if not grouping:
            raise core.Skip("No precomputed grouping to cleanup")
        return False

    def remove(self, *args, **kwargs):
        self.logger.info("Nothing to remove")

    def data(self, *args, **kwargs):
        self.cached('nr', remove=True)
        return None
