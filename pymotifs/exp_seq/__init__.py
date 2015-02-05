from pymotifs import core

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.mapping import Loader as MappingLoader
from pymotifs.exp_seq.positions import Loader as PositionsLoader


class Loader(core.MultiLoader):
    stages = [InfoLoader, PositionsLoader, MappingLoader]
