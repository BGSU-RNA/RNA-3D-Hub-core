from pymotifs import core

from pymotifs.autonomous.info import Loader as InfoLoader
from pymotifs.autonomous.chains import Loader as ChainLoader


class Loader(core.MultiStageLoader):
    stages = set([InfoLoader, ChainLoader])
