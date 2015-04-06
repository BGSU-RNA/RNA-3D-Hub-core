from pymotifs import core

from pymotifs.chains.info import Loader as InfoLoader
from pymotifs.chains.best import BestChainsAndModelsLoader
from pymotifs.chains.species import Loader as SpeciesLoader


class Loader(core.MultiStageLoader):
    dependencies = set([InfoLoader, BestChainsAndModelsLoader, SpeciesLoader])
