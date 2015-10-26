import pymotifs.core as core

from pymotifs.interactions.pairwise import Loader as PairwiseLoader
from pymotifs.interactions.helix_loop_summary import Loader as HLSummaryLoader


class Loader(core.StageContainer):
    stages = set([PairwiseLoader, HLSummaryLoader])
