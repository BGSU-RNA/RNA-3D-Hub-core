"""Load all PDB level infomation. This includes the folowing stages:

pdbs.info
    Fetch the basic information about a PDB.
"""

import pymotifs.core as core

from pymotifs.pdbs.info import Loader as InfoLoader


class Loader(core.StageContainer):
    """The loader for all PDB level data."""

    stages = set([InfoLoader])
    """Stages to run."""
