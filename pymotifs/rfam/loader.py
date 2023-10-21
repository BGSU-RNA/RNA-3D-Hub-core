"""Run all chain loading stages.
"""

from pymotifs import core

from pymotifs.rfam.map_to_rfam import Loader as MapLoader


class Loader(core.StageContainer):
    stages = set([MapLoader])
