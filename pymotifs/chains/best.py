from pymotifs import core
from pymotifs import models as mod
from pymotifs.mat_files import Loader as MatLoader


class BestChainsAndModelsLoader(core.SimpleLoader):
    dependencies = set([MatLoader])
    table = mod.PdbBestChainsAndModels

    def query(self, session, pdb, **kwargs):
        return session.query(mod.PdbBestChainsAndModels).filter_by(pdb_id=pdb)

    def data(self, pdb, **kwargs):
        # 'A,B,C', '1,2', ''
        mlab = core.Matlab(self.config['locations']['fr3d_root'])
        chains, models, err = mlab.loadBestChainsAndModels(pdb, nout=3)

        if err != '':
            print(chains, models, err)
            raise core.MatlabFailed(err)

        return {'pdb_id': pdb, 'best_chains': chains, 'best_models': models}
