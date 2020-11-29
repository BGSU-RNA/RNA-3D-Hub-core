"""Run all unit loaders. This will run the following unit stages:

- units.centers
- units.coordinates
- units.distances
- units.info
- units.quality
- units.rotation
"""

from pymotifs import core


# from pymotifs.md_units.rotation import Loader as MdRotationLoader

from pymotifs.md_units.centers import Loader as MdCentersLoader


class Loader(core.StageContainer):
    """A container for all unit level `Stage`s.
    """

    """The `Stages` this will run"""
    stages = set([MdCentersLoader])
