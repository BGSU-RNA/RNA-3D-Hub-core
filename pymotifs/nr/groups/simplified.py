import functools as ft
import itertools as it
import collections as coll

from pymotifs import core
from pymotifs.utils import result2dict
from pymotifs.constants import NR_DISCREPANCY_CUTOFF
from pymotifs.constants import EQUIVELANT_PAIRS

from pymotifs.models import PdbInfo
from pymotifs.models import ChainInfo
from pymotifs.models import ChainSpecies
from pymotifs.models import IfeInfo
from pymotifs.models import IfeChains
from pymotifs.models import ChainChainSimilarity
from pymotifs.models import ExpSeqInfo
from pymotifs.models import ExpSeqPdb

from pymotifs.utils import connectedsets as cs
from pymotifs.utils import correspondence as cr
from pymotifs.utils.structures import SYNTHEIC


def class_ordering(group):
    rep = group['members'][0]['chains'][0]
    return (rep['length'], rep['sequence'])


def ranking_key(chain):
    """Compute a key to order members of each group. The ordering produced
    should place the representative structure as the first one in the
    ordering.

    :chain: The chain to produce a key for.
    :returns: A value that can be used to sort the chains in descending
    order of quality.
    """

    bp_nt = 0.0
    if chain['bp']:
        bp_nt = float(chain['bp']) / float(chain['length'])

    name = list(ord(c) for c in chain['name'])
    return (-1 * bp_nt, -1 * chain['length'], name)


class Grouper(core.Base):
    """This is a class to determine the grouping of chains in a pdb. This will
    load all non-redundant chains and then group them by similarity. It
    produces a list of chains which are non-redundant. This does not process
    all chains as it only takes the best of each type from files. However, it
    does make a note about those other chains.
    """

    def valid_ife(self, ife):
        """Check if the given ife is valid. If not a warning statement will be
        produced and false returned. Ifes are valid if they have a non zero
        length. This constraint exists because some structures contain only
        modified RNA bases. These cannot be grouped (or really processed) by
        our nr set. So we remove them as to not get a series of 0 length groups
        in our nr clustering.

        :param dict ife: An ife as loaded by `self.ifes`.
        :returns: A flag to indicate if the ife is valid.
        """

        if ife['length']:
            return True
        self.logger.warning("Invalid ife %s, 0 length", ife['id'])
        return False

    def ifes(self, pdb):
        """Load all ife chains from a given pdb. This will get the RNA chains as
        well as load some interaction data about the chains.

        :pdb: A pdb to get the chains for.
        :returns: A list of dictionaries with data about all chains. The data
        is that which is provided by the info method.
        """

        with self.session() as session:
            query = session.query(ChainInfo.sequence,
                                  ChainInfo.chain_name.label('name'),
                                  IfeChains.chain_id.label('db_id'),
                                  IfeChains.is_integral,
                                  IfeChains.is_accompanying,
                                  IfeInfo.ife_id.label('id'),
                                  IfeInfo.bp_count.label('bp'),
                                  IfeInfo.pdb_id.label('pdb'),
                                  IfeInfo.length,
                                  PdbInfo.resolution,
                                  ChainSpecies.species_id.label('source')).\
                join(IfeInfo,
                     IfeInfo.pdb_id == ChainInfo.pdb_id).\
                join(IfeChains,
                     (IfeChains.ife_id == IfeInfo.ife_id) &
                     (IfeChains.chain_id == ChainInfo.chain_id)).\
                join(PdbInfo,
                     PdbInfo.pdb_id == ChainInfo.pdb_id).\
                join(ChainSpecies,
                     ChainSpecies.chain_id == ChainInfo.chain_id).\
                join(ExpSeqPdb,
                     ExpSeqPdb.chain_id == ChainInfo.chain_id).\
                join(ExpSeqInfo,
                     ExpSeqInfo.exp_seq_id == ExpSeqPdb.exp_seq_id).\
                filter(IfeInfo.pdb_id == pdb).\
                filter(IfeInfo.new_style == True).\
                filter(ExpSeqInfo.was_normalized == 1).\
                order_by(IfeChains.ife_id, IfeChains.index)

            if query.count() == 0:
                self.logger.warn("No ifes found for %s" % pdb)
                return []

            grouped = it.groupby(it.imap(result2dict, query),
                                 lambda g: g['id'])
            groups = []
            for ife_id, chains in grouped:
                chains = list(chains)
                groups.append({
                    'id':  ife_id,
                    'pdb': chains[0]['pdb'],
                    'bp': chains[0]['bp'],
                    'name': chains[0]['name'],
                    'length': chains[0]['length'],
                    'species': chains[0]['source'],
                    'chains': chains,
                    'db_id': chains[0]['db_id'],
                    'resolution': chains[0]['resolution'],
                })

        if not groups:
            raise core.InvalidState("No stored ifes for %s", pdb)

        self.logger.info("Found %i ifes for %s", len(groups), pdb)
        return groups

    def discrepancy(self, groups):
        chain_ids = []
        for group in groups:
            chain_ids.append(group['db_id'])

        with self.session() as session:
            query = session.query(ChainChainSimilarity).\
                filter(ChainChainSimilarity.chain_id_1.in_(chain_ids)).\
                filter(ChainChainSimilarity.chain_id_2.in_(chain_ids))

            discrepancy = coll.defaultdict(dict)
            for result in query:
                id1 = result.chain_id_1
                id2 = result.chain_id_2
                discrepancy[id1][id2] = result.discrepancy
                discrepancy[id2][id1] = result.discrepancy

        return dict(discrepancy)

    def alignments(self, chains):
        """Load the data about alignments between all pairs of chains.

        :chains: A list of chains.
        :returns: A dictionary of dictionaries.
        """

        pdbs = list(set(chain['pdb'] for chain in chains))
        helper = cr.Helper(self.config, self.session.maker)
        return helper.aligned_chains(pdbs)

    def has_good_alignment(self, all_alignments, group1, group2):
        """Check if the longest chain of the two groups align well.
        """

        db_id1 = group1['db_id']
        db_id2 = group2['db_id']
        if db_id1 not in all_alignments:
            self.logger.debug("Splitting %s %s: No alignments for %s",
                              group1['id'], group2['id'], group1['id'])
            return False

        alignments = all_alignments[db_id1]
        if db_id2 not in alignments:
            self.logger.debug("Splitting %s %s: No alignment between them",
                              group1['id'], group2['id'])
            return False

        if not alignments[db_id2]:
            self.logger.debug("Splitting %s %s: Bad alignment",
                              group1['id'], group2['id'])
            return False

        self.logger.debug("Good alignment: %s, %s", group1['id'], group2['id'])
        return True

    def are_similar_species(self, group1, group2):
        """Check if the longest chains of the two groups agree.
        """

        source1 = group1['source']
        source2 = group2['source']
        if source1 != SYNTHEIC[0] and source2 != SYNTHEIC[0] and \
                source1 is not None and source2 is not None \
                and source1 != source2:
            self.logger.debug("Splitting %s, %s: Different sources (%i, %i)",
                              group1['id'], group2['id'], source1, source2)
            return False

        self.logger.debug("Good species: %s, %s", group1['id'], group2['id'])
        return True

    def has_good_discrepancy(self, all_discrepancy, group1, group2):
        db_id1 = group1['db_id']
        db_id2 = group2['db_id']
        if not all_discrepancy:
            return True

        if db_id1 not in all_discrepancy:
            self.logger.warning("No computed discrepancy for %s", group1['id'])
            return True

        discrepancy = all_discrepancy[db_id1]
        if db_id2 not in discrepancy:
            self.logger.debug("Splitting %s %s: No discrepancy between them",
                              group1['id'], group2['id'])
            return False

        if discrepancy[db_id2] > NR_DISCREPANCY_CUTOFF:
            self.logger.debug("Splitting %s %s: Too large discrepancy %f",
                              group1['id'], group2['id'], discrepancy[db_id2])
            return False

        self.logger.debug("Good discrepancy: %s, %s",
                          group1['id'], group2['id'])
        return True

    def is_hard_coded_join(self, group1, group2):
        def simple(ife_id):
            parts = ife_id.split('|')
            return '|'.join([parts[0], parts[-1]])

        id1 = (simple(group1['id']), simple(group2['id']))
        id2 = (simple(group2['id']), simple(group1['id']))
        return id1 in EQUIVELANT_PAIRS or id2 in EQUIVELANT_PAIRS

    def are_equivalent(self, alignments, discrepancies, group1, group2):
        """Determine if two chains are equivalent. This requires that the
        chains align well and are not of conflicting species.

        :alignments: A dictionary of dictionaries about the alignment between
        chains.
        :chain1: The first chain.
        :chain2: The second chain.
        :returns: True if the two chains are equivalent, False otherwise
        """
        if self.is_hard_coded_join(group1, group2):
            return True

        return self.are_similar_species(group1, group2) and \
            self.has_good_alignment(alignments, group1, group2) and \
            self.has_good_discrepancy(discrepancies, group1, group2)

    def validate(self, connections, group):
        pairs = it.product(group, group)
        pairs = it.ifilter(lambda (a, b): a != b, pairs)
        pairs = it.ifilter(lambda (a, b): b not in connections[a], pairs)
        pairs = list(pairs)
        for pair in pairs:
            self.logger.debug("Pair %s, %s not connected", *pair)
        self.logger.debug("%i pairs are not connected", len(pairs))

    def pairs(self, chains, alignments, discrepancies):
        """Generate an iterator of all equivalent pairs of chains.

        :chains: The chains to build pairs of.
        :alignments: Alignments to use for checking validity.
        :discrepancies: Discrepancies of the chains to use.
        :returns: An iterable of all valid pairs.
        """

        mapping = {}
        for chain in chains:
            db_id = chain['db_id']
            if db_id in mapping:
                self.logger.error("db_id of entry already in mapping")
                self.logger.error("Entry: %s", chain)
                self.logger.error("Mapping: %s", mapping)
                raise core.InvalidState("Invalid mapping duplicated %i", db_id)
            mapping[db_id] = chain

        def is_mapped(index):
            def fn(p):
                known = p[index] in mapping
                if not known:
                    self.logger.warning("Unmapped db_id %i", p[index])
                return known
            return fn

        missing = it.imap(lambda c: c['db_id'], chains)
        missing = it.ifilterfalse(lambda c: c in alignments, missing)
        missing = it.product(missing, repeat=2)
        missing = list(missing)
        self.logger.info("Detected %i chains not in alignments", len(missing))

        equiv = ft.partial(self.are_equivalent, alignments, discrepancies)
        pairs = it.imap(lambda c: c['db_id'], chains)
        pairs = it.product(pairs, repeat=2)
        pairs = it.ifilter(is_mapped(0), pairs)
        pairs = it.ifilter(is_mapped(1), pairs)
        pairs = it.ifilter(lambda p: p[0] != p[1], pairs)
        pairs = it.imap(lambda p: (mapping[p[0]], mapping[p[1]]), pairs)
        pairs = it.ifilter(lambda p: p[0] == p[1] or equiv(*p), pairs)

        return pairs

    def group(self, chains, alignments, discrepancies):
        """Group all chains into connected components.

        :chains: A list of chains to group.
        :returns: A list of lists of the connected components.
        """

        mapping = {}
        graph = coll.defaultdict(set)
        for chain in chains:
            mapping[chain['id']] = chain
            graph[chain['id']].add(chain['id'])

        for chain1, chain2 in self.pairs(chains, alignments, discrepancies):
            self.logger.debug("Equivalent: %s %s", chain1['id'], chain2['id'])
            graph[chain1['id']].add(chain2['id'])
            graph[chain2['id']].add(chain1['id'])

        self.logger.info("Created %i connections", len(graph))
        groups = []
        for ids in cs.find_connected(graph).values():
            self.validate(graph, list(ids))
            groups.append(mapping[id] for id in ids)

        return groups

    def __call__(self, pdbs, **kwargs):
        """Group all chains in the given list of pdbs.

        :pdbs: A list of pdb ids to group the chains in.
        :returns: A list of groups dictionaries. These will contain one entry,
        'members' a list of chains in this group. These chains will be sorted
        by quality as defined by ranking_key. The chain data structure will
        be a dictionary as returned by merged_chains, but with an extra key
        'rank' which indicates the rank of the chain.
        """

        self.logger.info("Will group %i pdbs", len(pdbs))
        ifes = it.imap(self.ifes, pdbs)
        ifes = it.chain.from_iterable(ifes)
        ifes = it.ifilter(self.valid_ife, ifes)
        ifes = list(ifes)
        if not ifes:
            raise core.InvalidState("No ifes found in given pdbs")
        self.logger.info("Found %i ifes to cluster", len(ifes))

        alignments = self.alignments(ifes)
        if not alignments:
            raise core.InvalidState("No alignments loaded")

        discrepancy = self.discrepancy(ifes)
        if not discrepancy:
            self.logger.warning("No discrepancy data to cluster with")

        groups = []
        for group in self.group(ifes, alignments, discrepancy):
            members = []
            sorted_group = sorted(group, key=ranking_key)
            for index, member in enumerate(sorted_group):
                member['rank'] = index
                members.append(member)

            groups.append({'members': members})

        grouped_count = it.imap(lambda g: len(g['members']), groups)
        grouped_count = sum(grouped_count)
        if grouped_count != len(ifes):
            raise core.InvalidState("Only %i of %i ifes were grouped" %
                                    (grouped_count, len(ifes)))

        self.logger.info("Created %i groups with %i ifes",
                         len(groups), grouped_count)

        return groups
