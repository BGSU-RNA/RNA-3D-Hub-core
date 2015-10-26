"""This is the main module for running updates. Running this package will run
an update for all parts of the pipeline.
"""

from pymotifs import core
from pymotifs import download
from pymotifs import species_mapping
from pymotifs.units import loader as units
from pymotifs.pdbs import loader as pdbs
from pymotifs.loops import loader as loops
from pymotifs.chains import loader as chains
from pymotifs.interactions import loader as interactions
from pymotifs.exp_seq import loader as exp_seq
from pymotifs.correspondence import loader as correspondence
from pymotifs.ife import loader as ife
from pymotifs.export import loader as export


class Loader(core.StageContainer):
    stages = set([
        download.Downloader,
        pdbs.Loader,
        chains.Loader,
        species_mapping.Loader,
        units.Loader,
        interactions.Loader,
        loops.Loader,
        exp_seq.Loader,
        correspondence.Loader,
        ife.Loader,
        export.Loader
    ])
