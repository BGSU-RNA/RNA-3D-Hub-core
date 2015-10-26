from pymotifs import core

from pymotifs.loops.extractor import Loader as Extractor
from pymotifs.loops.positions import Loader as PositionLoader


class Loader(core.MultiStageLoader):
    stages = set([Extractor, PositionLoader])
