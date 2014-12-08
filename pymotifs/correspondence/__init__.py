import core

from correspondence.nts import Loader as NtLoader
# from correspondence.loops import Loader as LoopLoader
from correspondence.info import Loader as InfoLoader
from correspondence.interactions import Loader as InterLoader
# from correspondence.summary import Loader as SummaryLoader


class Loader(core.MultiLoader):
    stages = [InfoLoader, NtLoader, InterLoader]
    # stages = [InfoLoader, NtLoader, LoopLoader]
    # stages = [InfoLoader, NtLoader, LoopLoader, DiscLoader]
