"""
Download new files up to a date, process them completely, export, and chain_chain.comparison
After running this, it will be much faster to run the nr.loader stage and skip
all dependencies.
The intended use is to load years' worth of DNA data, then make weekly releases in a separate step
"""

from pymotifs import core
from pymotifs import download
from pymotifs.units import loader as units
from pymotifs.pdbs import loader as pdbs
from pymotifs.loops import loader as loops
from pymotifs.chains import loader as chains
from pymotifs.interactions import loader as interactions
from pymotifs.exp_seq import loader as exp_seq
from pymotifs.correspondence import loader as correspondence
from pymotifs.ife import loader as ife
from pymotifs.chain_chain import loader as chain_chain
from pymotifs.export import loader as export
from pymotifs.quality import loader as quality

class Loader(core.StageContainer):
    stages = set([
        download.Downloader,
        pdbs.Loader,
        chains.Loader,
        units.Loader,
        quality.Loader,
        interactions.Loader,
        loops.Loader,
        exp_seq.Loader,
        correspondence.Loader,
        ife.Loader,
        chain_chain.Loader,
        export.Exporter
    ])
