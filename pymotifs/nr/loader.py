"""Run all nr level stages. This runs:

nr.release
    Create a new NR release and cache the NR data.
nr.classes
    Store all NR classes (equivalence classes).
nr.parents
    Store the parent information for the cached NR release.
nr.parent_counts
    Store the counts of changes from the parents.
nr.id_mapping
    Modified the cached data with database ids.
nr.cqs
    Compute per-class Composite Quality Score data and
    store in the database.
nr.ordering
    Compute the ordering of chains within an nr class.
nr.cleanup
    Cleanup the cached data.
"""
import pymotifs.core as core

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.classes import Loader as ClassLoader
# from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.parents import Loader as ParentLoader
from pymotifs.nr.parent_counts import Loader as CountLoader
from pymotifs.nr.id_mapping import Loader as IdLoader
from pymotifs.nr.cqs import NrQualityLoader
from pymotifs.nr.ordering import Loader as OrderingLoader
from pymotifs.nr.cleanup import Cleanup
from pymotifs.nr.class_rank import Loader as ClassRankLoader

class Loader(core.StageContainer):
    # Temporarily skip the ordering stage because it is slow
    # stages = [ReleaseLoader, IdLoader, ClassLoader, ParentLoader,
    #           CountLoader, NrQualityLoader, Cleanup,ClassRankLoader]

    stages = [ReleaseLoader, IdLoader, ClassLoader, ParentLoader,
              CountLoader, NrQualityLoader, Cleanup,ClassRankLoader,OrderingLoader]
