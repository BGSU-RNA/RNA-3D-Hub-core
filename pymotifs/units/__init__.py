from pymotifs import core

from pymotifs.units.info import Loader as InfoLoader
# from pymotifs.units.coordinates import Loader as CoordinateLoader
# from pymotifs.units.distances import Loader as DistancesLoader
from pymotifs.units.quality import Loader as QualityLoader


class Loader(core.MultiLoader):
    stages = [InfoLoader, QualityLoader]
