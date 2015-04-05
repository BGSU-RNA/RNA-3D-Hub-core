from pymotifs import core
from pymotifs.models import ExpSeqInfo as Info
from pymotifs.models import ChainInfo
from pymotifs.chains.info import Loader as ChainLoader


class Loader(core.Loader):
    dependencies = set(ChainLoader)

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Info).\
                join(ChainInfo, ChainInfo.sequence == Info.sequence).\
                filter(ChainInfo.pdb_id == pdb)
            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Info.id).\
                join(ChainInfo, ChainInfo.sequence == Info.sequence).\
                filter(ChainInfo.pdb_id == pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            query = session.query(Info.id).\
                filter(Info.id.in_(ids)).\
                delete(synchronize_session=False)

    def sequences(self, pdb):
        with self.session() as session:
            query = session.query(ChainInfo.sequence).\
                filter_by(entity_macromolecule_type='Polyribonucleotide (RNA)').\
                filter_by(pdb_id=pdb)

            return set(result.sequence for result in query)

    def data(self, pdb, **kwargs):
        sequences = self.sequences(pdb)
        return [Info(sequence=seq, length=len(seq)) for seq in sequences]
