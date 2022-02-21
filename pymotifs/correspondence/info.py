"""Compute all new pairs of experimental sequences to compare. This will load
pairs for all given structures. It does not create comparisons to structures
not included with the given pdbs.
"""

import itertools as it
# import functools as ft
# from operator import itemgetter

from pymotifs import core
from pymotifs import utils as ut

from pymotifs.models import ExpSeqInfo
from pymotifs.models import ExpSeqPdb
from pymotifs.models import ChainSpecies
from pymotifs.models import CorrespondenceInfo

from pymotifs.constants import SYNTHENIC_SPECIES_ID
from pymotifs.constants import CORRESPONDENCE_HUGE_CUTOFF
from pymotifs.constants import CORRESPONDENCE_EXACT_CUTOFF

from pymotifs.exp_seq.loader import Loader as ExpSeqLoader
from pymotifs.chains.info import Loader as ChainInfoLoader
from pymotifs.chains.species import Loader as ChainSpeciesLoader


class Loader(core.MassLoader):
    """A class to load up all pairs of experimental sequences that should have
    an alignment attempted. This does not work per structure as many other
    things do, but instead will compute all possible pairs and then work per
    pair, inserting or storing each pair as needed.
    """

    dependencies = set([ChainSpeciesLoader, ExpSeqLoader, ChainInfoLoader])
    allow_no_data = True
    table = CorrespondenceInfo

    exact_cutoff = CORRESPONDENCE_EXACT_CUTOFF                                      ## = 36 """Length cutoff before being matched as small"""
    huge_cutoff = CORRESPONDENCE_HUGE_CUTOFF                                        ## = 2000 """Length cutoff before being matched as huge"""

    def has_data(self, pdb, **kwargs):
        """Always returns false. This stage is odd in that we do the filtering
        later on, so we always skip to computing data.
        """
        return False

    def lookup_sequences(self, pdb):
        """Return all exp_seq_ids for the given pdb. This only assign the
        species id from the given pdb.

        :param str pdb: The pdb id to get all sequences for.
        :returns: A list of dictionaries of unique sequences.
        """

        with self.session() as session:
            query = session.query(ExpSeqPdb.exp_seq_id.label('id'),                 ## There are some label functions, they are here for rename the column to "id", "length", and "species" 
                                  ExpSeqInfo.normalized_length.label('length'),
                                  ChainSpecies.species_id.label('species'),
                                  ExpSeqInfo.entity_type.label('entity_type')).\
                join(ExpSeqInfo,
                     ExpSeqInfo.exp_seq_id == ExpSeqPdb.exp_seq_id).\
                outerjoin(ChainSpecies,
                          ChainSpecies.chain_id == ExpSeqPdb.chain_id).\
                filter(ExpSeqPdb.pdb_id == pdb).\
                filter(ExpSeqInfo.was_normalized).\
                distinct()

            if not query.count():
                self.logger.warning("No sequences for %s" % pdb)
            
            for result in query:
                self.logger.info("show the query result:%s %s %s %s" % (result.id, result.length, result.species, result.entity_type))
                self.logger.info(type(result.species))
                self.logger.info("show the lookup_sequences return: %s" % ut.row2dict(result))                  ## the return value example: {'length': 76L, 'id': 4177L, 'species': 4932L}
                
            return [ut.row2dict(result) for result in query]

    def length_match(self, pair):
        """Determine if the length of the sequences in the pair match.

        Parameters
        ----------
        pair : (dict, dict)
            The pair of sequences to compare.

        Returns
        -------
        length_matches : bool
            If the lengths match
        """

        smallest = min(p['length'] for p in pair)
        largest = max(p['length'] for p in pair)

        if smallest < self.exact_cutoff:
            return smallest == largest

        lower = max(0.5 * largest, self.exact_cutoff)
        upper = 2 * smallest

        if smallest >= self.huge_cutoff:
            lower = self.huge_cutoff
        else:
            upper = min(upper, self.huge_cutoff)

        return lower <= smallest <= largest <= upper

    def species_matches(self, pair):
        """Check if a pair of sequences have compatible species. This means the
        same species or having None and synthetic (32360) assigned.

        :param tuple pair: The pair of sequences to compare.
        :returns: The
        """

        common = pair[0]['species'].intersection(pair[1]['species'])
        self.logger.info("show the current species: %s  %s" % (pair[0]['entity_type'], pair[1]['entity_type']))
        if len(common):
            return True
        species = set()
        for seq in pair:
            species.update(seq['species'])
        return len(species) <= 1 or None in species or SYNTHENIC_SPECIES_ID in species

    @property
    def known(self):
        """Get all known pairs. This will look up the known pairs the first
        time it is used, otherwise it will return the same set of known pairs.

        :returns: A set of all known pairs. The entries in the set will be of
        the form (exp_seq_id_1, exp_seq_id_2).
        """

        if not hasattr(self, '_known'):
            with self.session() as session:
                query = session.query(
                    self.table.exp_seq_id_1.label('id1'),
                    self.table.exp_seq_id_2.label('id2'),
                ).distinct()

                self._known = set((result.id1, result.id2) for result in query)
        return self._known

    def is_known(self, pair):
        """Determine if this pair has already been computed.

        :param tuple pair: The pair of sequences.
        :returns: A boolean.
        """

        return (pair[0]['id'], pair[1]['id']) in self.known

    def ids_with_column_names(self, pair):
        """Turn a tuple into a dictionary for the database.

        :param tuple pair: The sequence pair.
        :returns: A dictionary.
        """

        return {
            'exp_seq_id_1': pair[0]['id'],
            'exp_seq_id_2': pair[1]['id']
        }

    def unique_sequences(self, sequences):
        """Merge all sequences into unique sequences. It will merge all species
        into a single set.

        :param list sequences: A list of the sequences loaded from
        `lookup_sequences`.
        :returns: A list of unique sequences.
        """

        mapping = {}
        for seq in sequences:
            sid = seq['id']
            self.logger.info("show the seq info: %s %s %s %s" % (seq["id"], seq["length"],seq["species"],seq["entity_type"]))
            if sid not in mapping:
                mapping[sid] = dict(seq)
                mapping[sid]['species'] = set([seq['species']])
            mapping[sid]['species'].update([seq['species']])
            self.logger.info("show the mapping info %s" % mapping)
        return mapping.values()

    def entity_type_matches(self, pair):
        """Check if a pair of entity_types have compatible species. 

        :param tuple pair: The pair of entity_types to compare.
        :returns: The boolean
        """
        if pair[0]['entity_type'] == 'dna':
            self.logger.info("the entity_types of the current pair %s %s" % (pair[0]['entity_type'], pair[1]['entity_type']))
        return pair[0]['entity_type'] == pair[1]['entity_type']


    def is_match(self, pair):
        """Check if a pair of sequences are a pair. A pair is a match if they
        are a good length, have matching sequences and is not known.

        :param tuple pair: A pair of sequences to compare.
        :returns: A boolean.
        """

        return self.length_match(pair) and \
            self.species_matches(pair) and \
            self.entity_type_matches(pair) and \
            not self.is_known(pair)

    def sequences(self, pdbs):
        """Get all unique_sequences for all given pdbs.

        :param list pdbs: The pdbs to lookup.
        :returns: A list of unique sequences.
        """

        self.logger.info("Using %i pdbs", len(pdbs))
        seqs = it.imap(self.lookup_sequences, pdbs)
        self.logger.info("show the it.map's result: %s" % seqs)
        seqs = it.chain.from_iterable(seqs)
        seqs = self.unique_sequences(seqs)
        self.logger.info("show the unique_sequences's result: %s" % seqs)
        if not seqs:
            raise core.InvalidState("Found no new sequences")

        self.logger.info("Found %i unique sequences", len(seqs))
        return sorted(seqs, key=lambda s: s['id'])

    def pairs(self, pdbs):
        """Compute all pairs for the given pdbs

        :param list pdbs: The list pdbs to process.
        :returns: The list pairs
        """

        seqs = self.sequences(pdbs)                                         ## seqs is a list of dicts with information about sequences.                       
        self.logger.info("show the seqs info of pairs function: %s" % seqs[0]) ## [{'species': set(4932), 'length': 76, 'id': 4177, 'entity_type': 'rna'}]
        # pairs = it.combinations(seqs, 2)                                    ## just a guess, the function will return a list of diff combinations if we have more than 2 seqs. Thus, if we only have one seq, it will return [].
        # Thus, the problem is here, we can add condiction for each entity type.
        rna_pairs = it.combinations([f for f in seqs if f['entity_type']=='rna'],2)
        dna_pairs = it.combinations([f for f in seqs if f['entity_type']=='dna'],2)
        hybrid_pairs = it.combinations([f for f in seqs if f['entity_type']=='hybrid'],2)
        # pairs = list(rna_pairs)+list(dna_pairs)+list(hybrid_pairs)
        # self.logger.info("Found pairs %s", pairs)                           ## <itertools.combinations object at 0x7fb8f641dec0>
        self_pairs = it.izip(seqs, seqs)                                    ## Ex. ({'species': {4932}, 'length': 76, 'id': 4177, 'entity_type': 'rna'}, {'species': {4932}, 'length': 76, 'id': 4177, 'entity_type': 'rna'})
        self.logger.info("Found self_pairs %s", self_pairs)                 ## <itertools.izip object at 0x7fb8f6405bd8>
        # return sorted(it.chain.from_iterable([self_pairs, pairs]),
        #               key=lambda p: (p[0]['id'], p[1]['id']))
        
        return sorted(it.chain.from_iterable([self_pairs,rna_pairs,dna_pairs,hybrid_pairs]),
                      key=lambda p: (p[0]['id'], p[1]['id']))        

    def data(self, pdbs, **kwargs):
        """Compute all new correspondences pairs.

        :param list pdbs: The pdbs to process.
        :returns: A list of pairs to store.
        """

        pairs = self.pairs(pdbs)
        pairs = it.ifilter(self.is_match, pairs)
        pairs = sorted(pairs, key=lambda p: (p[0]['id'], p[1]['id']))
        pairs = it.imap(self.ids_with_column_names, pairs)
        pairs = list(pairs)
        self.logger.info("Found %i new correspondence pairs", len(pairs))
        return pairs
