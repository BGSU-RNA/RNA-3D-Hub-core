"""
This is the main module for running updates. Running this package will run an
update for all parts of the pipeline.
"""

from pymotifs import core
from pymotifs import download
from pymotifs import units
from pymotifs import pdbs
# from pymotifs import loops
from pymotifs import chains
from pymotifs import interactions
from pymotifs import exp_seq
from pymotifs import correspondence
from pymotifs import species_mapping
from pymotifs import ife


class Loader(core.MultiStageLoader):
    dependencies = set([
        download.Downloader,
        pdbs.Loader,
        chains.Loader,
        species_mapping.Loader,
        units.Loader,
        interactions.Loader,
        # loops.Loader,
        exp_seq.Loader,
        correspondence.Loader,
        ife.Loader,
    ])
