from pymotifs import core

from pymotifs.models import ExpSeqInfo as Exp
from pymotifs.models import ExpSeqPosition as Position


class Loader(core.MassLoader):
    insert_max = 5000

    def sequences(self):
        with self.session() as session:
            query = session.query(Exp).\
                outerjoin(Position, Position.exp_seq_id == Exp.id).\
                filter(Position.id == None)
            return [(result.id, result.sequence) for result in query]

    def positions(self, exp_id, sequence):
        positions = []
        for index, char in enumerate(sequence):
            positions.append({'exp_seq_id': exp_id,
                              'unit': char,
                              'index': index})
        return positions

    def data(self, pdbs, **kwargs):
        for (exp_id, sequence) in self.sequences():
            for position in self.positions(exp_id, sequence):
                yield Position(**position)
