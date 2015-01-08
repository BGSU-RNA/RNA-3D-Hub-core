import pymotifs.core as core

from pymotifs.interactions.pairwise import Loader as PairwiseLoader
# from pymotifs.interactions.summary import Loader as SummaryLoader


class InteractionLoader(core.MultiLoader):
    stages = [PairwiseLoader]
