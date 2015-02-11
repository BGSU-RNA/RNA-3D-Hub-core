from pymotifs import core

from pymotifs.chains.info import Loader as InfoLoader
from pymotifs.chains.best import BestChainsAndModelsLoader


class Loader(core.MultiStageLoader):
    stages = [InfoLoader, BestChainsAndModelsLoader]
