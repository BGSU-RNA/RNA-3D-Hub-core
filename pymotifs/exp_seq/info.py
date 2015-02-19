from pymotifs import core
from pymotifs.models import ExpSeqInfo as Info
from pymotifs.models import ChainInfo


class Loader(core.MassLoader):

    def known(self):
        with self.session() as session:
            query = session.query(Info.sequence)
            return set(result.sequence for result in query)

    def possible(self):
        with self.session() as session:
            query = session.query(ChainInfo.sequence).\
                filter_by(entity_macromolecule_type='Polyribonucleotide (RNA)')

            return set(result.sequence for result in query)

    def data(self, *args, **kwargs):
        possible = self.possible()
        known = self.known()
        new = possible - known
        return [Info(sequence=seq, length=len(seq)) for seq in new]
