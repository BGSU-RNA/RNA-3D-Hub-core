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
        """Fetch the (sequence,type) pairs to process . This will use the given pdbs to
        extract all sequences come from nucleic acid chains in those structures. This
        will then get all the unique sequences paired with the macromolecule type.

        Parameters
        ----------
        pdbs : list
            The list of pdb ids to use.

        Returns
        -------
        seq_type_pairs : list
            A list of (sequence,type) pairs to process.
            For example, ('AGAACAUUC','rna')
        """

        self.logger.info("to_process method")

        simplify_type = {}
        simplify_type['Polydeoxyribonucleotide (DNA)'] = 'dna'
        simplify_type['Polyribonucleotide (RNA)'] = 'rna'
        simplify_type['polyribonucleotide'] = 'rna'
        simplify_type['polydeoxyribonucleotide'] = 'dna'
        simplify_type['polydeoxyribonucleotide/polyribonucleotide hybrid'] = 'hybrid'
        simplify_type['DNA/RNA Hybrid'] = 'hybrid'

        macromolecule_types = simplify_type.keys()

        data = set([])   # set so we get unique (sequence,type) pairs

        with self.session() as session:

            query = session.query(mod.ChainInfo.sequence,
                                  mod.ChainInfo.entity_macromolecule_type).\
                filter(mod.ChainInfo.pdb_id.in_(pdbs)).\
                filter(mod.ChainInfo.entity_macromolecule_type.in_(macromolecule_types))   
            data.update((result.sequence,simplify_type[result.entity_macromolecule_type]) for result in query)
        # sort the set of pairs into a list, first by sequence length, then by sequence, then by type
        return sorted(data, key=lambda p: (len(p[0]), p[0], p[1]))  # return a list of pairs

    def identify_hybrid_sequences(self):
        """
            This function is trying to get all sequences whose entity_macromolecule_type is hybrid
        """
        """
        Returns
        -------
            A list of sequences
        """
        hybrid_sequences = []
        with self.session() as session:
            query = session.query(mod.ChainInfo.sequence).\
                                  filter(mod.ChainInfo.entity_macromolecule_type.in_(['polydeoxyribonucleotide/polyribonucleotide hybrid',\
                                  'DNA/RNA Hybrid'])).distinct()
            for i in query:
                hybrid_sequences.append(i[0])

        return(hybrid_sequences)

    def query(self, session, seq_type):
        """The query to find and remove all exp seq info entries for a given
        sequence.

        Parameters
        ----------
        session : pymotifs.core.session.Session
            The session object to use.
        seq_type : tuple
            The (sequence,entity_type) to look up.
        """

        sequence = seq_type[0]
        entity_type = seq_type[1]

        return session.query(mod.ExpSeqInfo).filter_by(md5=self.md5(sequence)).filter_by(entity_type=entity_type)

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

        # There may be a more sophisticated translation method available, then use that
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

        # New on 10/13/2021, treat T as an acceptable character, store rna/dna/hybrid
        if character in set(['A', 'C', 'G', 'T', 'U', 'N']):
            return character
        return self.translation.get(character, None)

    def normalize(self, sequence, entity_type):
        """Normalize a sequence. This will translate all characters in the
        sequence so that ...

        Parameters
        ----------
        sequence : str
            The sequence to normalize
        entity_type : str
            rna, dna, or hybrid

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

    def data(self, seq_type, **kwargs):
        """Compute data to store for the given sequence. This will compute what
        is needed to store, including the md5 hash which should be unique
        across all sequences.

        Parameters
        ----------
        seq_type = (seq,type) : pair of strings
            The sequence to store, RNA/DNA/hybrid entity type

        Returns
        -------
        data : dict
            A dictonary with 'sequence', 'md5', 'normalized', 'length',
            'normalized_length', 'was_normalized', 'entity_type' entries.
        """

        # split pair into sequence and type
        seq = seq_type[0]
        entity_type = seq_type[1]

        self.logger.info("Sequence is %s" % seq)
        self.logger.info("Entity type is %s" % entity_type)

        # convert unusual letters
        normalized = self.normalize(seq, entity_type)


        return {
            'sequence': seq,
            'md5': self.md5(seq),
            'normalized': normalized,
            'length': len(seq),
            'normalized_length': len(normalized) if normalized else 0,
            'was_normalized': normalized is not None,
            'entity_type': entity_type
        }
