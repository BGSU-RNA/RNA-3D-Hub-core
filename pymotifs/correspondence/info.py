from pymotifs import core
from pymotifs import utils as ut

from pymotifs.models import ExpSeqInfo
from pymotifs.models import ExpSeqPdb
from pymotifs.models import ChainSpecies
from pymotifs.models import CorrespondenceInfo as Info
from pymotifs.models import ExpSeqChainMapping

from pymotifs.exp_seq.info import Loader as ExpSeqInfoLoader
from pymotifs.exp_seq.chain_mapping import Loader as ExpSeqChainMappingLoader
from pymotifs.chains.info import Loader as ChainInfoLoader
from pymotifs.chains.species import Loader as ChainSpeciesLoader


class Loader(core.Loader):
    """A class to load up all pairs of experimental sequences that should have
    an alignment attempted. This does not work per structure as many other
    things do, but instead will compute all possible pairs and then work per
    pair, inserting or storing each pair as needed.
    """

    dependencies = set([ChainSpeciesLoader, ExpSeqInfoLoader,
                        ExpSeqChainMappingLoader, ChainInfoLoader])

    allow_no_data = True
    mark = False

    def has_data(self, pdb, **kwargs):
        return False

    def remove(self, pdb, **kwargs):
        self.logger.info("We never automatically remove correspondence info")
        pass

    def lookup_sequences(self, pdb):
        """Return all exp_seq_ids for the given pdb.
        """
        with self.session() as session:
            query = session.query(ExpSeqPdb.exp_seq_id,
                                  ExpSeqInfo.length,
                                  ChainSpecies.species_id).\
                join(ExpSeqInfo,
                     ExpSeqInfo.exp_seq_id == ExpSeqPdb.exp_seq_id).\
                outerjoin(ChainSpecies,
                          ChainSpecies.chain_id == ExpSeqPdb.chain_id).\
                filter(ExpSeqPdb.pdb_id == pdb)

            return [ut.row2dict(result) for result in query]

    def pairs(self, exp_seq):
        """Compute the pairs of sequence ids which should be aligned. This does
        not check if those pairs have already been aligned, it just computes
        ones to align. This will find pairs of sequences which have similar
        length and are from organisms which do not conflict. That means it will
        align things which have no known species or are synthetic.
        """

        id1 = exp_seq['exp_seq_id']
        length = exp_seq['length']
        species = exp_seq['species_id']

        with self.session() as session:
            query = session.query(ExpSeqInfo).\
                filter(ExpSeqInfo.exp_seq_id != id1)

            # We treat large and small sequences differently, for small
            # sequences (< 36 nts) we have to have an exact match. For large
            # sequences we require that the sequences be of similar length,
            # which is from 0.5 to two times the first. This is a broad enough
            # range to cover all good alignments.
            if length < 36:
                query = query.filter(ExpSeqInfo.length == length)
            else:
                shortest = max(0.5 * length, 36)
                length_terms = ((ExpSeqInfo.length >= shortest) &
                                (ExpSeqInfo.length <= 2 * length))

                if length >= 2000:
                    length_terms &= (ExpSeqInfo.length >= 2000)
                else:
                    length_terms &= (ExpSeqInfo.length < 2000)
                query = query.filter(length_terms)

            # We also only get pairs between things that have
            # 32360 is the taxon id for synthetic
            if species is not None and species != 32360:
                query = query.\
                    join(ExpSeqChainMapping,
                         ExpSeqChainMapping.exp_seq_id == ExpSeqInfo.exp_seq_id).\
                    join(ChainSpecies,
                         ChainSpecies.chain_id == ExpSeqChainMapping.chain_id).\
                    filter((ChainSpecies.species_id.in_([32360, species])) |
                           (ChainSpecies.species_id == None))

            pairs = [{'exp_seq_id1': id1, 'exp_seq_id2': id1}]
            for result in query.distinct():
                id2 = result.exp_seq_id
                pairs.append({
                    'exp_seq_id1': min(id1, id2),
                    'exp_seq_id2': max(id1, id2)
                })

            return sorted(pairs,
                          key=lambda e: (e['exp_seq_id1'], e['exp_seq_id2']))

    def is_known(self, pair):
        with self.session() as session:
            ids = [pair['exp_seq_id1'], pair['exp_seq_id2']]
            query = session.query(Info).\
                filter(Info.exp_seq_id_1.in_(ids)).\
                filter(Info.exp_seq_id_2.in_(ids))
            return bool(query.limit(1).count())

    def to_info(self, pair):
        return Info(exp_seq_id_1=pair['exp_seq_id1'],
                    exp_seq_id_2=pair['exp_seq_id2'])

    def data(self, pdb, **kwargs):
        exp_seqs = self.lookup_sequences(pdb)
        data = []
        for exp_seq in exp_seqs:
            pairs = self.pairs(exp_seq)
            data.extend(self.to_info(p) for p in pairs if not self.is_known(p))

        if not data:
            raise core.Skip("No possible pairings found for %s" % pdb)
        return data
