"""Determine best and chains and models for structures.

Runs matlab across the given pdbs files to determine which are the best chains
and models. These are used for the motif atlas to determine which chains and
models to use for extracting motifs.
"""

from pymotifs import core
from pymotifs.utils import matlab
from pymotifs import models as mod
from pymotifs.mat_files import Loader as MatLoader


class BestChainsAndModelsLoader(core.SimpleLoader):
    dependencies = set([MatLoader])
    @property
    def table(self):
        return mod.PdbBestChainsAndModels

    def query(self, session, pdb, **kwargs):
        return session.query(mod.PdbBestChainsAndModels).filter_by(pdb_id=pdb)

    def data(self, pdb, **kwargs):
        # 'A,B,C', '1,2', ''
        mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
        chains, models, err = mlab.loadBestChainsAndModels(pdb, nout=3)

        if err != '':
            raise matlab.MatlabFailed(err)

        return {'pdb_id': pdb, 'best_chains': chains, 'best_models': models}
