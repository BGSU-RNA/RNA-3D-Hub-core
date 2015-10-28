import hashlib

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import grouper
from pymotifs.utils.structures import Structure
from pymotifs.chains.info import Loader as ChainLoader


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader])
    mark = False

    def to_process(self, pdbs, **kwargs):
        sequences = set()
        helper = Structure(self.session.maker)
        for chunk in grouper(100, pdbs):
            chain_ids = helper.rna_chains(pdbs, return_id=True)
            chain_ids = set(p[1] for p in chain_ids)
            with self.session() as session:
                query = session.query(mod.ChainInfo.sequence).\
                    filter(mod.ChainInfo.chain_id.in_(chain_ids)).\
                    distinct()

                for result in query:
                    sequences.add(result.sequence)

        return list(sequences)

    def query(self, session, sequence):
        return session.query(mod.ExpSeqInfo).\
            filter(mod.ExpSeqInfo.md5 == self.md5(sequence))

    def md5(self, sequence):
        return hashlib.md5(sequence).hexdigest()

    def data(self, seq, **kwargs):
        return mod.ExpSeqInfo(sequence=seq, md5=self.md5(seq), length=len(seq))
