"""Run all quality loaders. This will run the following quality stages:

- quality.downloader
- quality.units
- quality.pdb
"""

from pymotifs import core

from pymotifs.quality.download import Loader as Downloader
from pymotifs.quality.units import Loader as UnitsQuality
from pymotifs.quality.pdb import Loader as PdbQuality
from pymotifs.quality.clashes import Loader as ClashLoader


class Loader(core.StageContainer):
    """A container for all unit level `Stage`s.
    """

    """The `Stages` this will run"""
    stages = set([Downloader, UnitsQuality, PdbQuality, ClashLoader])
