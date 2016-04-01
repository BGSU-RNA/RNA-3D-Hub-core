from pymotifs import core

from pymotifs.units.info import Loader as InfoLoader
# from pymotifs.units.coordinates import Loader as CoordinateLoader
from pymotifs.units.distances import Loader as DistancesLoader
from pymotifs.units.quality import Loader as QualityLoader
from pymotifs.units.redundant import RedundantNucleotidesLoader
from pymotifs.units.centers import Loader as CenterLoader


class Loader(core.StageContainer):
    stages = set([InfoLoader, QualityLoader, DistancesLoader, CenterLoader
                  RedundantNucleotidesLoader])
