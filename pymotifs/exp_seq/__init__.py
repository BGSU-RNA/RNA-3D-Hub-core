import core

from exp_seq.info import Loader as InfoLoader
from exp_seq.mapping import Loader as MappingLoader
from exp_seq.positions import Loader as PositionsLoader


class Loader(core.MultiLoader):
    stages = [InfoLoader, PositionsLoader, MappingLoader]
