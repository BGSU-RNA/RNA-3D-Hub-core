"""
Run all unit loaders. This will run the following unit stages:

- units.centers
- units.coordinates
- units.distances
- units.info
- units.quality
- units.rotation
"""

from pymotifs import core

from pymotifs.units.info import Loader as InfoLoader
from pymotifs.units.coordinates import Loader as CoordinateLoader
from pymotifs.units.center_rotation import Loader as CenterRotationLoader
from pymotifs.units.incomplete import Loader as IncompleteLoader
from pymotifs.units.bond_orientation import Loader as BondOrientationLoader


class Loader(core.StageContainer):
    """
    A container for all unit level stages
    """

    """
    The stages this will run
    """
    stages = set([InfoLoader, CenterRotationLoader,
                  CoordinateLoader, IncompleteLoader, BondOrientationLoader])
