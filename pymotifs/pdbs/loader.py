"""
Load all PDB level infomation. This includes the folowing stages:

pdbs.info:    Fetch the basic information about a PDB.
assemblies:   Which chains and symmetries are in which assemblies
"""

import pymotifs.core as core

from pymotifs.pdbs.info import Loader as InfoLoader
from pymotifs.pdbs.assemblies import Loader as AssemblyLoader


class Loader(core.StageContainer):
    """The loader for all PDB level data."""

    stages = set([InfoLoader,AssemblyLoader])
