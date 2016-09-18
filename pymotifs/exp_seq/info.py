"""Load data to the exp_seq_info table. This will process the given PDBs and
extract the unique experimental sequences from then and write data about each
sequence to the exp_seq_info table.
"""

import hashlib

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import grouper
from pymotifs.utils.structures import Structure
from pymotifs.chains.info import Loader as ChainLoader


class Loader(core.SimpleLoader):
    """The actual loader for this stage."""

    dependencies = set([ChainLoader])
    mark = False

    @property
    def table(self):
        return mod.ExpSeqInfo

    def to_process(self, pdbs, **kwargs):
        """Fetch the sequences to process . This will use the given pdbs to
        extract all sequences come from RNA chains in those structures. This
        will then get all the unique sequences.

        Parameters
        ----------
        pdbs : list
            The list of pdb ids to use.

        Returns
        -------
        sequences : list
            A list of sequences to process.
        """

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

        Parameters
        ----------
        session : pymotifs.core.session.Session
            The session object to use.
        sequence : str
            The sequence to lookup.
        """
        return session.query(mod.ExpSeqInfo).filter_by(md5=self.md5(sequence))

    def md5(self, sequence):
        """Compute the md5 hash of a sequence.

        Parameters
        ----------
        sequence : str
            The sequence to get the MD5 of

        Returns
        -------
        md5 : str
            The md5 hash as a hex digest.
        """
        return hashlib.md5(sequence).hexdigest()

    @property
    def translation(self):
        """The translation dictonary for this loader. It will translate
        characters from modified notation to A, C, G, U only. It will lookup
        the entries in rna_unit_modified_correspondencies as well as harded
        codes: {'X': 'N', 'I': 'G'}.

        Returns
        -------
        translation : dict
            A dictonary that maps from non standard units to standard units.
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
        """Translate sequences to a standard representation. If no standard
        representation is known then it will return None.

        Parameters
        ----------
        character : str
            A character to translate to A, C, G, U or N.

        Returns
        -------
        translated : str
            The translated character.
        """
        if character in set(['A', 'C', 'G', 'U', 'N']):
            return character
        return self.translation.get(character, None)

    def normalize(self, sequence):
        """Normalize a sequence. This will translate all characters in the
        sequence so that

        Parameters
        ----------
        sequence : str
            The sequence to normalize

        Returns
        -------
        normalized : str
            The normalized sequence.
        """

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
        """Compute data to store for the given sequence. This will compute what
        is needed to store, including the md5 hash which should be unique
        across all sequences.

        Parameters
        ----------
        seq : str
            The sequence to store.

        Returns
        -------
        data : dict
            A dictonary with 'sequence', 'md5', 'normalized', 'length',
            'normalized_length' and 'was_normalized' entries.
        """

        normalized = self.normalize(seq)
        return {
            'sequence': seq,
            'md5': self.md5(seq),
            'normalized': normalized,
            'length': len(seq),
            'normalized_length': len(normalized) if normalized else 0,
            'was_normalized': normalized is not None
        }
