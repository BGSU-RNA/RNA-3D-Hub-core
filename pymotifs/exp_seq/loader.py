"""Load all experimental sequence data. This will run:

exp_seq.info
    Get all information about experimental sequences
exp_seq.mapping
    Map units to all experimental sequence positions
exp_seq.positions
    Compute all positions in the experimental sequences
exp_seq.chain_mapping
    Map chains to experimental sequences
"""

from pymotifs import core

from pymotifs.exp_seq.info import Loader as InfoLoader
from pymotifs.exp_seq.mapping import Loader as MappingLoader
from pymotifs.exp_seq.positions import Loader as PositionsLoader
from pymotifs.exp_seq.chain_mapping import Loader as ChainMappingLoader


class Loader(core.StageContainer):
    stages = set([InfoLoader, ChainMappingLoader, PositionsLoader,
                  MappingLoader])
