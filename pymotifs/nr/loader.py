import pymotifs.core as core

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.classes import Loader as ClassLoader
from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.parents import Loader as ParentLoader
from pymotifs.nr.parent_counts import Loader as CountLoader
from pymotifs.nr.id_mapping import Loader as IdLoader
from pymotifs.nr.cleanup import Cleanup


class Loader(core.StageContainer):
    stages = [ReleaseLoader, IdLoader, ClassLoader, ChainLoader, ParentLoader,
              CountLoader, Cleanup]
