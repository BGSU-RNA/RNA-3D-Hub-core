import pymotifs.core as core

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.classes import Loader as ClassLoader
from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.parents import Loader as ParentLoader


class Loader(core.StageContainer):
    stages = [ReleaseLoader, ClassLoader, ChainLoader, ParentLoader]
