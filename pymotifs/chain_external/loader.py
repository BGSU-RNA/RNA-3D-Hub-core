"""Run all chain_external level stages:

- mapping
"""

from pymotifs import core
from pymotifs.chain_external.mapping import Loader as chain_external_mapping_loader


class Loader(core.StageContainer):
    """Contains all chain external level stages
    """

    """The stages to run."""
    stages = set([chain_external_mapping_loader])
