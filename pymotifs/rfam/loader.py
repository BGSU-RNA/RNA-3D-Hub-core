"""
Run all stages that map chains to Rfam
Order:
map_to_rfam
consensus_name
"""

from pymotifs import core

from pymotifs.chains.taxid_species_domain import Loader as DomainLoader
from pymotifs.rfam.map_to_rfam import Loader as MapLoader      # map PDB chains to Rfam families, align them
from pymotifs.rfam.consensus_name import Loader as NameLoader  # assign domain, Rfam, consensus name to PDB chains


class Loader(core.StageContainer):
    stages = set([MapLoader])
    stages = set([NameLoader])
    stages = set([DomainLoader,MapLoader,NameLoader])
