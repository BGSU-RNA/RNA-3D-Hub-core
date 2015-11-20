import itertools as it

from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod

from pymotifs.exp_seq.info import Loader as ExpSeqInfoLoader
from pymotifs.exp_seq.chain_mapping import Loader as ExpSeqChainMappingLoader
from pymotifs.chains.info import Loader as ChainInfoLoader
from pymotifs.chains.species import Loader as ChainSpeciesLoader


class Loader(core.SimpleLoader):
    """A class to load up all pairs of experimental sequences that should have
    an alignment attempted. This does not work per structure as many other
    things do, but instead will compute all possible pairs and then work per
    pair, inserting or storing each pair as needed.
    """

    dependencies = set([ChainSpeciesLoader, ExpSeqInfoLoader,
                        ExpSeqChainMappingLoader, ChainInfoLoader])

    mark = False
    table = mod.CorrespondenceInfo

    small_cutoff = 36
    huge_cutoff = 2000

    def valid_length(self, pair):
        """We treat large and small sequences differently, for small
        sequences (< 36 nts) we have to have an exact match. For large
        sequences we require that the sequences be of similar length,
        which is from 0.5 to two times the first. This is a broad enough
        range to cover all good alignments.
        """

        length1 = min(pair[0]['length'], pair[1]['length'])
        length2 = max(pair[0]['length'], pair[1]['length'])
        if length1 < self.small_cutoff:
            return length1 == length2

        lower_limit = max(int(0.5 * length2), self.small_cutoff)
        if length2 > self.huge_cutoff or length1 > self.huge_cutoff:
            lower_limit = self.huge_cutoff

        upper_limit = 2 * length1
        if length1 < self.huge_cutoff or length2 < self.huge_cutoff:
            upper_limit = 2000

        return lower_limit <= length1 <= length2 <= upper_limit

    def valid_species(self, pair):
        """Check weather a pair has valid sequences. A pair has valid sequences
        if the species match or the species do not conflict. This happens when
        the species is either None or 32360, the species id for synthetic
        constructs.
        """
        species = it.imap(lambda p: p['species'], pair)
        species = set(it.chain.from_iterable(species))
        return None in species or 32360 in species or len(species) == 1

    def pairs(self, entries):
        """Produce all possible pairs of sequences from the given entries. This
        will create unique pairs where the species do not conflict and the
        length's are close enough to produce a good alignment. The result will
        be sorted by the ids so the order will always be the same.
        """

        pairs = it.product(entries, repeat=2)
        pairs = it.ifilter(lambda p: p[0]['id'] <= p[1]['id'], pairs)
        pairs = it.ifilter(self.valid_length, pairs)
        pairs = it.ifilter(self.valid_species, pairs)
        return sorted(pairs, key=lambda p: (p[0]['id'], p[1]['id']))

    def exp_seqs(self, pdbs):
        """Look up all experimental sequences for
        """
        with self.session() as session:
            query = session.query(mod.ExpSeqPdb.exp_seq_id.label('id'),
                                  mod.ExpSeqInfo.length,
                                  mod.ChainSpecies.species_id.label('species')).\
                join(mod.ExpSeqInfo,
                     mod.ExpSeqInfo.exp_seq_id == mod.ExpSeqPdb.exp_seq_id).\
                join(mod.ExpSeqChainMapping,
                     mod.ExpSeqChainMapping.exp_seq_id == mod.ExpSeqPdb.exp_seq_id).\
                join(mod.ChainSpecies,
                     mod.ChainSpecies.chain_id == mod.ExpSeqChainMapping.chain_id).\
                filter(mod.ExpSeqPdb.pdb_id.in_(pdbs)).\
                distinct().\
                order_by(mod.ExpSeqInfo.exp_seq_id)

            return [ut.row2dict(result) for result in query]

    def merge(self, exps):
        merged = []
        key = lambda d: d['id']
        for exp_id, entries in it.groupby(sorted(exps, key=key), key):
            entries = list(entries)
            merged.append({
                'id': exp_id,
                'length': entries[0]['length'],
                'species': set(e['species'] for e in entries)
            })
        return merged

    def to_process(self, pdbs, **kwargs):
        """Compute the pairs of experimental sequences to process. This will
        look up all possible pairs in the given list of pdbs and then filter
        them to only those which match our length and species cutoffs.
        """
        exps = self.exp_seqs(pdbs)
        merged = self.merge(exps)
        return self.pairs(merged)

    def query(self, session, pair):
        return session.query(self.table).\
            filter_by(exp_seq_id_1=pair[0]['id'], exp_seq_id_2=pair[1]['id'])

    def data(self, pair, **kwargs):
        return {
            'exp_seq_id_1': pair[0]['id'],
            'exp_seq_id_2': pair[1]['id'],
        }
