"""
Run the motif atlas loaders.
We typically run them separately, loop mapping every week, motif atlas every 4 weeks
"""

from pymotifs import core

from pymotifs.motif_atlas.mapping import Loader as MappingLoader
from pymotifs.motif_atlas.make_atlas import Loader as MakeAtlasLoader


class Loader(core.StageContainer):
    """
    A container for all unit level stages
    """

    """
    The stages this will run
    """
    stages = set([MappingLoader, MakeAtlasLoader])
