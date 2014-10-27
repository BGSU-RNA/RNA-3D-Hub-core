import logging

import core
import utils

from models import ExpSeqInfo as Exp
from models import ExpSeqPosition as Position

from rnastructure.tertiary.cif import CIF

logger = logging.getLogger(__name__)


class Loader(core.Loader):
    name = 'exp_seq_positions'
    update_gap = False
    insert_max = 5000
    allow_no_data = True

    def __init__(self, config, maker):
        self.finder = utils.CifFileFinder(config)
        super(Loader, self).__init__(config, maker)

    def has_data(self, pdb):
        with self.session() as session:
            query = self.__query__(session, pdb)
            return bool(query.count())

    def remove(self, pdb):
        with self.session() as session:
            query = session.query(Exp.id).filter(Exp.pdb == pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            session.query(Position).\
                filter(Position.exp_seq_id.in_(ids)).\
                delete(synchronize_session=False)

    def transform(self, pdb):
        mapped = []
        with self.session() as session:
            query = session.query(Exp).filter_by(pdb=pdb)
            for result in query:
                mapped.append((result.id, result.pdb, result.chain))
        return mapped

    def data(self, entry, **kwargs):
        exp_id, pdb, chain = entry
        with open(self.finder(pdb), 'rb') as raw:
            cif = CIF(raw)

        data = []
        seen = set()
        chain = cif.chain('1_555', 1, chain)
        seq = chain.experimental_sequence_mapping()
        for index, (seq, seq_id, _) in enumerate(seq):
            if seq_id not in seen:
                seen.add(seq_id)
                data.append(Position(id=seq_id, unit=seq,
                                     index=index + 1,
                                     exp_seq_id=exp_id))
        return data

    def __query__(self, session, pdb):
        return session.query(Position).\
            join(Exp, Exp.id == Position.exp_seq_id).\
            filter(Exp.pdb == pdb)
