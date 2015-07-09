from pymotifs import core
from pymotifs.chain_chain.comparision import Loader as CompareLoader


class Loader(core.MultiStageLoader):
    dependencies = set([CompareLoader])
