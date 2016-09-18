"""Run all interaction level stages. This runs:

interactions.pairwise
    Annotate all pairwise interactions
interactions.flanking
    Annotate all flanking interactions
interactions.summary
    Summarize the number of interactions for each unit.
"""

import pymotifs.core as core

from pymotifs.interactions.pairwise import Loader as PairwiseLoader
from pymotifs.interactions.flanking import Loader as FlankingLoader
from pymotifs.interactions.summary import Loader as SummaryLoader


class Loader(core.StageContainer):
    stages = set([PairwiseLoader, FlankingLoader, SummaryLoader])
