from pymotifs import core
from pymotifs import utils

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.models import ExpSeqPosition as Position

from rnastructure.tertiary.cif import CIF


class Loader(core.Loader):
    name = 'exp_seq_positions'
    update_gap = False
    insert_max = 5000
    allow_no_data = True

    def __init__(self, config, maker):
        self.finder = utils.CifFileFinder(config)
        super(Loader, self).__init__(config, maker)

    def has_data(self, entry):
        exp_id, pdb, chain = entry
        with self.session() as session:
            query = session.query(Position).filter_by(exp_seq_id=exp_id)
            return bool(query.count())

    def remove(self, entry):
        with self.session() as session:
            session.query(Position).filter_by(exp_seq_id=entry[0]).\
                delete(synchronize_session=False)

    def transform(self, pdb, **kwargs):
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
        try:
            seq = chain.experimental_sequence_mapping()
        except ValueError as err:
            self.logger.warning("Can't map %s, %s", pdb, chain)
            raise core.SkipValue(str(err))

        for index, (seq, seq_id, _) in enumerate(seq):
            if seq_id not in seen:
                seen.add(seq_id)
                data.append(Position(id=seq_id, unit=seq,
                                     index=index + 1,
                                     exp_seq_id=exp_id))
        return data
