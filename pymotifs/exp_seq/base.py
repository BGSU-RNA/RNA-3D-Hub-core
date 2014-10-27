import core
import utils
from models import ObsSequence as Obs

from rnastructure.tertiary.cif import CIF


class Base(core.Loader):
    klass = None

    def __init__(self, config, maker):
        self.finder = utils.CifFileFinder(config)
        super(Base, self).__init__(config, maker)

    def remove(self, pdb):
        with self.session() as session:
            query = session.query(Exp.id).filter(Exp.pdb == pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            session.query(klass).\
                filter(klass.exp_seq_id.in_(ids)).\
                delete(synchronize_session=False)

    def transform(self, pdb):
        mapped = []
        with self.session() as session:
            query = session.query(Exp).filter_by(pdb=pdb)
            for result in query:
                mapped.append((result.id, result.pdb, result.chain))
        return mapped
