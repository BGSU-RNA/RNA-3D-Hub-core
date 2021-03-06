"""Run all chain chain level stages:

- comparison
"""

from pymotifs import core
from pymotifs.chain_chain.comparison import Loader as CompareLoader


class Loader(core.StageContainer):
    """Contains all chain chain level stages
    """

    """The stages to run."""
    stages = set([CompareLoader])
