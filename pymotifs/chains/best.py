from pymotifs import core
from pymotifs.models import PdbBestChainsAndModels
from pymotifs.matfiles import Loader as MatLoader


class BestChainsAndModelsLoader(core.SimpleLoader):
    dependencies = set(MatLoader)

    def __init__(self, *args, **kwargs):
        super(BestChainsAndModelsLoader, self).__init__(*args, **kwargs)
        self.matlab = core.Matlab(self.config['locations']['fr3d_root'])

    def query(self, session, pdb, **kwargs):
        return session.query(PdbBestChainsAndModels). \
            filter(PdbBestChainsAndModels.pdb_id == pdb)

    def data(self, pdb, **kwargs):
        # 'A,B,C', '1,2', ''
        best_chains, best_models, err_msg = \
            self.matlab.loadBestChainsAndModels(str(pdb), nout=3)

        if err_msg != '':
            raise core.MatlabFailed(err_msg)

        return PdbBestChainsAndModels(pdb_id=pdb, best_chains=best_chains,
                                      best_models=best_models)
