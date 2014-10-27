import core
import utils
from models import ExpSeqInfo as Info

from rnastructure.tertiary.cif import CIF


class Loader(core.Loader):
    name = 'exp_seq_info'
    update_gap = False

    def __init__(self, config, maker):
        self.finder = utils.CifFileFinder(config)
        super(Loader, self).__init__(config, maker)

    def remove(self, pdb):
        with self.session() as session:
            session.query(Info).filter(Info.pdb == pdb).delete()

    def has_data(self, pdb):
        with self.session() as session:
            count = session.query(Info).filter(Info.pdb == pdb).count()
        return bool(count)

    def data(self, pdb, **kwargs):
        data = []
        with open(self.finder(pdb), 'rb') as raw:
            cif = CIF(raw)
            for chain in cif.chains():
                # TODO: Sequence should be single letter form
                # sequence = ''.join(chain.experimental_sequence())
                data.append(Info(pdb=pdb, chain=chain['chain']))
        return data
