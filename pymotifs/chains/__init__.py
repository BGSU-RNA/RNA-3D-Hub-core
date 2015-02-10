from pymotifs import core

from pymotifs.chains.info import Loader as InfoLoader


class Loader(core.MultiStageLoader):
    stages = [InfoLoader]
