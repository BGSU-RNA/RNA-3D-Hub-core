from pymotifs import core
from pymotifs.models import ExpSeqInfo as Info


class Loader(core.SimpleLoader):
    name = 'exp_seq_info'
    update_gap = False

    def __init__(self, config, maker):
        super(Loader, self).__init__(config, maker)

    def query(self, session, pdb):
        return session.query(Info).filter(Info.pdb == pdb)

    def data(self, pdb, **kwargs):
        seen = set()
        data = []
        for chain in self.structure(pdb).chains():
            if (pdb, chain['chain']) not in seen:
                seen.add((pdb, chain['chain']))
                data.append(Info(pdb=pdb, chain=chain['chain']))
        return data
