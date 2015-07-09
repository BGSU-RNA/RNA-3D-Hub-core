import pymotifs.core as core

from pymotifs.pdbs.info import Loader as InfoLoader
from pymotifs.pdbs.obsolete import Loader as ObsoleteLoader


class Loader(core.MultiStageLoader):
    dependencies = set([InfoLoader, ObsoleteLoader])
