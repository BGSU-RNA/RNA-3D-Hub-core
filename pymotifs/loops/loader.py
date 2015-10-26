from pymotifs import core

from pymotifs.loops.extractor import Loader as Extractor
from pymotifs.loops.positions import Loader as PositionLoader


class Loader(core.StageContainer):
    stages = set([Extractor, PositionLoader])
