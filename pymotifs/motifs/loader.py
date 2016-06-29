import pymotifs.core as core

from pymotifs.motifs.cleanup import Loader as Cleanup
from pymotifs.motifs.motifs import Loader as MotifLoader
# from pymotifs.motifs.parents import Loader as ParentLoader
from pymotifs.motifs.release import Loader as ReleaseLoader
from pymotifs.motifs.loop_order import Loader as LoopOrderLoader
from pymotifs.motifs.discrepancies import Loader as DiscrepancyLoader
from pymotifs.motifs.loop_positions import Loader as LoopPositionLoader
# from pymotifs.motifs.parent_counts import Loader as ParentCountsLoader


class Loader(core.StageContainer):
    stages = [ReleaseLoader, DiscrepancyLoader, MotifLoader, LoopOrderLoader,
              LoopPositionLoader, Cleanup]
