from pymotifs import core
from pymotifs.loops.extractor import Loader as Extractor


class Loader(core.MultiStageLoader):
    stages = [Extractor]
