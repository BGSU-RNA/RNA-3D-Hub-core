from pymotifs import core

from pymotifs.correspondence.nts import Loader as NtLoader
# from correspondence.loops import Loader as LoopLoader
from pymotifs.correspondence.info import Loader as InfoLoader
from pymotifs.correspondence.interactions import Loader as InterLoader
# from correspondence.summary import Loader as SummaryLoader


class Loader(core.MultiLoader):
    stages = [InfoLoader, NtLoader, InterLoader]
    # stages = [InfoLoader, NtLoader, LoopLoader]
    # stages = [InfoLoader, NtLoader, LoopLoader, DiscLoader]
