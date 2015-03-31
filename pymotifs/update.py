"""
This is the main module for running updates. Running this package will run an
update for all parts of the pipeline.
"""

from pymotifs import core
from pymotifs import download
from pymotifs import export
from pymotifs import units
from pymotifs import pdbs
from pymotifs import loops
from pymotifs import chains
from pymotifs import interactions


class Loader(core.MultiStageLoader):
    stages = [download.Downloader,
              export.CifAtom,
              pdbs.Loader,
              chains.Loader,
              units.Loader,
              interactions.Loader,
              loops.Loader]
