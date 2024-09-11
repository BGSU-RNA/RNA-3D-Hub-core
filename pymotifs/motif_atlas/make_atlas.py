"""Load all motif data.

This will create a new release, if needed, and then store all data.
"""

import pymotifs.core as core

from pymotifs.motif_atlas.cleanup import Loader as Cleanup
from pymotifs.motif_atlas.info import Loader as MotifLoader
from pymotifs.motif_atlas.parents import Loader as ParentLoader
from pymotifs.motif_atlas.release import Loader as ReleaseLoader
from pymotifs.motif_atlas.loop_order import Loader as LoopOrderLoader
from pymotifs.motif_atlas.assignments import Loader as AssignmentLoader
from pymotifs.motif_atlas.discrepancies import Loader as DiscrepancyLoader
from pymotifs.motif_atlas.loop_positions import Loader as LoopPositionLoader
# from pymotifs.motif_atlas.parent_counts import Loader as ParentCountsLoader
from pymotifs.motif_atlas.signature import Loader as SignatureLoader  # bp_signature
#from pymotifs.motif_atlas.annotations import Loader as AnnotationLoader  # old annotations
from pymotifs.motif_atlas.secondary_structure import Loader as SecondaryStructureLoader


class Loader(core.StageContainer):
    stages = set([ReleaseLoader, DiscrepancyLoader, MotifLoader, LoopOrderLoader,
              LoopPositionLoader, AssignmentLoader, SignatureLoader,
              SecondaryStructureLoader, ParentLoader])

