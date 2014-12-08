import core
import utils
from models import ExpSeqInfo as Info

from rnastructure.tertiary.cif import ComplexOperatorException


class Loader(core.Loader):
    name = 'exp_seq_info'
    update_gap = False

    def __init__(self, config, maker):
        self.finder = utils.CifData(config)
        super(Loader, self).__init__(config, maker)

    def remove(self, pdb):
        with self.session() as session:
            session.query(Info).filter(Info.pdb == pdb).delete()

    def has_data(self, pdb):
        with self.session() as session:
            count = session.query(Info).filter(Info.pdb == pdb).count()
        return bool(count)

    def data(self, pdb, **kwargs):
        try:
            cif = self.finder(pdb)
        except utils.MissingFileException:
            raise core.SkipValue("Missing cif file %s" % pdb)
        except ComplexOperatorException:
            raise core.SkipValue("Skipping %s with complex operator" % pdb)

        seen = set()
        data = []
        for chain in cif.chains():
            if (pdb, chain['chain']) not in seen:
                seen.add((pdb, chain['chain']))
                data.append(Info(pdb=pdb, chain=chain['chain']))
        return data
