import itertools as it
import functools as ft
from operator import itemgetter

from pymotifs import core
from pymotifs import utils as ut

from pymotifs.models import ExpSeqInfo
from pymotifs.models import ExpSeqPdb
from pymotifs.models import ChainSpecies
from pymotifs.models import CorrespondenceInfo

from pymotifs.constants import CORRESPONDENCE_HUGE_CUTOFF
from pymotifs.constants import CORRESPONDENCE_SMALL_CUTOFF

from pymotifs.chains.info import Loader as ChainInfoLoader
from pymotifs.exp_seq.info import Loader as ExpSeqInfoLoader
from pymotifs.chains.species import Loader as ChainSpeciesLoader
from pymotifs.exp_seq.chain_mapping import Loader as ExpSeqChainMappingLoader


class Loader(core.MassLoader):
    """A class to load up all pairs of experimental sequences that should have
    an alignment attempted. This does not work per structure as many other
    things do, but instead will compute all possible pairs and then work per
    pair, inserting or storing each pair as needed.
    """

    dependencies = set([ChainSpeciesLoader, ExpSeqInfoLoader,
                        ExpSeqChainMappingLoader, ChainInfoLoader])
    allow_no_data = True
    table = CorrespondenceInfo

    small_cutoff = CORRESPONDENCE_SMALL_CUTOFF
    huge_cutoff = CORRESPONDENCE_HUGE_CUTOFF

    def has_data(self, pdb, **kwargs):
        return False

    def lookup_sequences(self, pdb):
        """Return all exp_seq_ids for the given pdb.
        """

        with self.session() as session:
            query = session.query(ExpSeqPdb.exp_seq_id.label('id'),
                                  ExpSeqInfo.normalized_length.label('length'),
                                  ChainSpecies.species_id.label('species')).\
                join(ExpSeqInfo,
                     ExpSeqInfo.exp_seq_id == ExpSeqPdb.exp_seq_id).\
                outerjoin(ChainSpecies,
                          ChainSpecies.chain_id == ExpSeqPdb.chain_id).\
                filter(ExpSeqPdb.pdb_id == pdb).\
                filter(ExpSeqInfo.was_normalized).\
                distinct()

            return [ut.row2dict(result) for result in query]

    def length_match(self, pair):
        smallest = min(p['length'] for p in pair)
        largest = max(p['length'] for p in pair)

        if smallest < self.small_cutoff:
            return smallest == largest

        lower = max(0.5 * largest, self.small_cutoff)
        upper = 2 * smallest

        if smallest >= self.huge_cutoff:
            lower = self.huge_cutoff
        else:
            upper = min(upper, self.huge_cutoff)

        return lower <= smallest <= largest <= upper

    def species_matches(self, seqs):
        species = set(s['species'] for s in seqs)
        return len(species) <= 1 or None in species or 32630 in species

    def is_known(self, known, pair):
        return (pair[0]['id'], pair[1]['id']) in known

    def known(self):
        with self.session() as session:
            query = session.query(
                self.table.exp_seq_id_1.label('id'),
                self.table.exp_seq_id_2.label('id'),
            ).distinct()

            return set((result.id1, result.id2) for result in query)

    def as_pair(self, pair):
        return {
            'exp_seq_id_1': pair[0]['id'],
            'exp_seq_id_2': pair[1]['id']
        }

    def unique_sequences(self, sequences):
        mapping = dict((seq['id'], seq) for seq in sequences)
        return mapping.values()

    def data(self, pdbs, **kwargs):
        self.logger.info("Using %i pdbs", len(pdbs))
        seqs = it.imap(self.lookup_sequences, pdbs)
        seqs = it.chain.from_iterable(seqs)
        seqs = sorted(self.unique_sequences(seqs), key=itemgetter('id'))
        self.logger.info("Found %i unique sequences", len(seqs))

        is_known = ft.partial(self.is_known, self.known())
        pairs = it.combinations(seqs, 2)
        pairs = it.ifilter(self.length_match, pairs)
        pairs = it.ifilter(self.species_matches, pairs)
        pairs = it.ifilterfalse(is_known, pairs)
        pairs = it.imap(self.as_pair, pairs)
        pairs = list(pairs)

        self.logger.info("Found %i new correspondence pairs", len(pairs))
        return pairs
