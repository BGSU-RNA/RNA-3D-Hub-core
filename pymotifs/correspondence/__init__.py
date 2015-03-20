from pymotifs import core

from pymotifs.correspondence.nts import Loader as NtLoader
from pymotifs.correspondence.info import Loader as InfoLoader
from pymotifs.correspondence.interactions import Loader as InterLoader


class Loader(core.MultiLoader):
    stages = [InfoLoader, NtLoader, InterLoader]
