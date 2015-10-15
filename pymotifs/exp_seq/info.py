import hashlib

from pymotifs import core
from pymotifs.utils import grouper
from pymotifs.utils.structures import Structure
from pymotifs.models import ExpSeqInfo as Info
from pymotifs.models import ChainInfo
from pymotifs.chains.info import Loader as ChainLoader


class Loader(core.Loader):
    dependencies = set([ChainLoader])

    mark = False

    def to_process(self, pdbs, **kwargs):
        sequences = set()
        helper = Structure(self.session.maker)
        for chunk in grouper(100, pdbs):
            chain_ids = helper.rna_chains(pdbs, return_id=True)
            chain_ids = set(p[1] for p in chain_ids)
            with self.session() as session:
                query = session.query(ChainInfo.sequence).\
                    filter(ChainInfo.chain_id.in_(chain_ids)).\
                    distinct()

                for result in query:
                    sequences.add(result.sequence)

        return list(sequences)

    def remove(self, sequence, **kwargs):
        md5 = self.md5(sequence)
        with self.session() as session:
            session.query(Info).\
                filter(Info.md5 == md5).\
                delete(synchronize_session=False)

    def has_data(self, sequence, **kwargs):
        md5 = self.md5(sequence)
        with self.session() as session:
            query = session.query(Info).\
                filter(Info.md5 == md5)
            return bool(query.count())

    def md5(self, sequence):
        return hashlib.md5(sequence).hexdigest()

    def data(self, seq, **kwargs):
        return Info(sequence=seq, md5=self.md5(seq), length=len(seq))
