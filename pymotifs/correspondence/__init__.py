import core

from correspondence.nts import Loader as NtLoader
from correspondence.loops import Loader as LoopLoader
from correspondence.info import Loader as InfoLoader
# from correspondence.discrepancy import Loader as DiscLoader


class Loader(core.MultiLoader):
    stages = [InfoLoader, NtLoader, LoopLoader]
    # stages = [InfoLoader, NtLoader, LoopLoader, DiscLoader]
