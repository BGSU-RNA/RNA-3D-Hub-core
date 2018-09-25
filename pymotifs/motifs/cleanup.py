"""Clean up all cached and release data.
"""

import os

import pymotifs.core as core

from pymotifs.motifs.info import Loader as MotifLoader
from pymotifs.motifs.release import Loader as ReleaseLoader
from pymotifs.motifs.loop_order import Loader as LoopOrderLoader
from pymotifs.motifs.discrepancies import Loader as DiscrepancyLoader
from pymotifs.motifs.loop_positions import Loader as LoopPositionLoader


class Loader(core.Loader):
    dependencies = set([MotifLoader, ReleaseLoader, LoopOrderLoader,
                        DiscrepancyLoader, LoopPositionLoader])
    allow_no_data = True

    def has_data(self, *args, **kwargs):
        for loop_type in ReleaseLoader.loop_types:
            if os.path.exists(self.cache_filename(loop_type)):
                return True

    def remove(self, *args, **kwargs):
        self.logger.info("Nothing to remove for motifs cleanup")

    def data(self, *arges, **kwargs):
        for loop_type in ReleaseLoader.loop_types:
            self.evict(loop_type)
