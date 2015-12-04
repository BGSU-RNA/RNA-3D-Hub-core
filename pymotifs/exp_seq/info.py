import hashlib

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import grouper
from pymotifs.utils.structures import Structure
from pymotifs.chains.info import Loader as ChainLoader


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader])
    mark = False
    table = mod.ExpSeqInfo

    def to_process(self, pdbs, **kwargs):
        sequences = set()
        helper = Structure(self.session.maker)
        for chunk in grouper(1000, pdbs):
            ids = set(p[1] for p in helper.rna_chains(chunk, return_id=True))
            with self.session() as session:
                query = session.query(mod.ChainInfo.sequence).\
                    filter(mod.ChainInfo.chain_id.in_(ids)).\
                    distinct()

                sequences.update(result.sequence for result in query)

        return sorted(sequences, key=lambda s: (len(s), s))

    def query(self, session, sequence):
        """The query to find and remove all exp seq info entries for a given
        sequence.
        """
        return session.query(mod.ExpSeqInfo).filter_by(md5=self.md5(sequence))

    def md5(self, sequence):
        """Compute the md5 hash of a sequence.
        """
        return hashlib.md5(sequence).hexdigest()

    @property
    def translation(self):
        """The translation dictonary for this loader. It will translate
        characters from modified notation to A, C, G, U only.
        """

        if hasattr(self, '_translation'):
            return self._translation

        # X is really N, but old PDB data doesn't respect that.
        self._translation = {'X': 'N', 'I': 'G'}
        table = mod.RnaUnitModifiedCorrespondencies
        with self.session() as session:
            query = session.query(
                table.rna_unit_modified_correspondencies_id,
                table.standard_unit,
            )

            for result in query:
                self._translation[result[0]] = result[1]
        return self._translation

    def translate(self, character):
        if character in set(['A', 'C', 'G', 'U', 'N']):
            return character
        return self.translation.get(character, None)

    def normalize(self, sequence):
        normalized = []
        size = len(sequence) - 1
        for index, seq in enumerate(sequence):
            norm = self.translate(seq)
            if not norm:
                if index == size:
                    self.logger.warning("Skipping final unit %s tRNA/AA", norm)
                    continue
                return None
            normalized.append(norm)

        return ''.join(normalized)

    def data(self, seq, **kwargs):
        normalized = self.normalize(seq)
        return {
            'sequence': seq,
            'md5': self.md5(seq),
            'normalized': normalized,
            'length': len(seq),
            'normalized_length': len(normalized) if normalized else 0,
            'was_normalized': normalized is not None
        }
