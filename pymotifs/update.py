"""Run the update pipeline.

This stage represents the update pipeline. Running this will run all stages
that are part of the pipeline proper.
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
from pymotifs.chain_chain import loader as cc
from pymotifs.export import loader as export
from pymotifs.nr import loader as nr
from pymotifs.quality import loader as quality


class Loader(core.StageContainer):
    stages = set([
        download.Downloader,
        pdbs.Loader,
        chains.Loader,
        species_mapping.Loader,
        units.Loader,
        quality.Loader,
        interactions.Loader,
        loops.Loader,
        exp_seq.Loader,
        correspondence.Loader,
        ife.Loader,
        cc.Loader,
        export.Exporter,
        nr.Loader,
    ])
