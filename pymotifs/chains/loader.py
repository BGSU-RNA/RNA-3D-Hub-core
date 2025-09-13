"""Run all chain loading stages.
"""

from pymotifs import core

from pymotifs.chains.info import Loader as InfoLoader
from pymotifs.chains.taxid_species_domain import Loader as TaxLoader

class Loader(core.StageContainer):
    stages = set([InfoLoader,TaxLoader])
