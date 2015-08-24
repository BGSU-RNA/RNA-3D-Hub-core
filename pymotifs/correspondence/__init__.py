from pymotifs import core

from pymotifs.correspondence.positions import Loader as PositionLoader
from pymotifs.correspondence.info import Loader as InfoLoader
from pymotifs.correspondence.cleanup import Loader as Cleanup
from pymotifs.correspondence.summary import Loader as SummaryLoader
from pymotifs.exp_seq import Loader as ExpSeqLoader
# from pymotifs.correspondence.interactions import Loader as InterLoader


class Loader(core.MultiStageLoader):
    stages = set([ExpSeqLoader, InfoLoader, PositionLoader, SummaryLoader,
                  Cleanup])
