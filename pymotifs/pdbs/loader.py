"""Load all PDB level infomation. This includes the folowing stages:

pdbs.info
    Fetch the basic information about a PDB.
pdbs.obsolete
    Fetch the data on obsolete pdbs.
"""

import pymotifs.core as core

from pymotifs.pdbs.info import Loader as InfoLoader
from pymotifs.pdbs.obsolete import Loader as ObsoleteLoader


class Loader(core.StageContainer):
    """The loader for all PDB level data."""

    stages = set([InfoLoader, ObsoleteLoader])
    """Stages to run."""
