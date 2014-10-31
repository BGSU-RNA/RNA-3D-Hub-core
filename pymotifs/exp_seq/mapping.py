import logging

import core
import utils

from models import ExpSeqInfo as Exp
from models import ExpSeqUnitMapping as Mapping

from rnastructure.tertiary.cif import CIF

logger = logging.getLogger(__name__)


class Loader(core.Loader):
    name = 'exp_seq_mapping'
    update_gap = False
    insert_max = 5000
    allow_no_data = True

    def __init__(self, config, maker):
        self.finder = utils.CifFileFinder(config)
        super(Loader, self).__init__(config, maker)

    def has_data(self, pdb):
        with self.session() as session:
            query = session.query(Mapping).\
                join(Exp, Exp.id == Mapping.exp_seq_id).\
                filter(Exp.pdb == pdb)
            return bool(query.count())

    def remove(self, pdb):
        with self.session() as session:
            query = session.query(Exp.id).filter(Exp.pdb == pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            session.query(Mapping).\
                filter(Mapping.exp_seq_id.in_(ids)).\
                delete(synchronize_session=False)

    def transform(self, pdb):
        mapped = []
        with self.session() as session:
            query = session.query(Exp).filter_by(pdb=pdb)
            for result in query:
                mapped.append((result.id, result.pdb, result.chain))
        return mapped

    def data(self, entry, **kwargs):
        obs_id, pdb, chain = entry
        mapping = []
        seen = set()
        with open(self.finder(pdb), 'rb') as raw:
            cif = CIF(raw)

        chain = cif.chain('1_555', 1, chain)
        if not chain:
            logging.error("Chain %s is unmappable", chain)
            return []

        seq = chain.experimental_sequence_mapping()
        for (_, seq_id, unit_id) in seq:
            if unit_id not in seen:
                mapping.append(Mapping(unit_id=unit_id,
                                       exp_seq_position_id=seq_id,
                                       exp_seq_id=obs_id))
                seen.add(unit_id)

        return mapping
