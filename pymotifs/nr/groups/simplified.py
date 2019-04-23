"""This is a module to compute a simple NR grouping. This grouping method uses
the basic approach where each chains must have similar sequence, similar source
organism and similar geometry. Other groups are possible, but this is simple
one that shows the behavior we generally like.

For details of this approach read:
"""

import functools as ft
import itertools as it
import collections as coll
import operator as op

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import result2dict
from pymotifs.utils import discrepancy as disc
from pymotifs.constants import NR_DISCREPANCY_CUTOFF
from pymotifs.constants import EQUIVALENT_PAIRS
from pymotifs.constants import SYNTHENIC_SPECIES_ID
from pymotifs.constants import NR_MIN_HOMOGENEOUS_SIZE

from pymotifs.utils import connectedsets as cs
from pymotifs.utils import correspondence as cr
from pymotifs.utils.structures import SYNTHEIC


def ranking_key(chain):
    """Compute a key to order members of each group. The ordering produced
    should place the representative structure as the first one in the
    ordering. This is the initial ordering which may be changed later when
    computing the representative.

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

    """A flag to determine if we should use discrepancies when splitting"""
    use_discrepancy = True

    """A flag to determine if we should use sepcies when splitting"""
    use_species = True

    """A flag to force groups to have distinct species"""
    must_enforce_single_species = True

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
            query = session.query(mod.ChainInfo.sequence,
                                  mod.ChainInfo.chain_name.label('name'),
                                  mod.IfeChains.chain_id.label('db_id'),
                                  mod.IfeChains.is_integral,
                                  mod.IfeChains.is_accompanying,
                                  mod.IfeInfo.ife_id.label('id'),
                                  mod.IfeInfo.bp_count.label('bp'),
                                  mod.IfeInfo.pdb_id.label('pdb'),
                                  mod.IfeInfo.length,
                                  mod.PdbInfo.resolution,
                                  mod.PdbInfo.experimental_technique.label('method'),
                                  mod.ChainSpecies.species_id.label('species')).\
                join(mod.IfeInfo,
                     mod.IfeInfo.pdb_id == mod.ChainInfo.pdb_id).\
                join(mod.IfeChains,
                     (mod.IfeChains.ife_id == mod.IfeInfo.ife_id) &
                     (mod.IfeChains.chain_id == mod.ChainInfo.chain_id)).\
                join(mod.PdbInfo,
                     mod.PdbInfo.pdb_id == mod.ChainInfo.pdb_id).\
                join(mod.ChainSpecies,
                     mod.ChainSpecies.chain_id == mod.ChainInfo.chain_id).\
                join(mod.ExpSeqPdb,
                     mod.ExpSeqPdb.chain_id == mod.ChainInfo.chain_id).\
                join(mod.ExpSeqInfo,
                     mod.ExpSeqInfo.exp_seq_id == mod.ExpSeqPdb.exp_seq_id).\
                filter(mod.IfeInfo.pdb_id == pdb).\
                filter(mod.IfeInfo.new_style == True).\
                filter(mod.ExpSeqInfo.was_normalized == 1).\
                order_by(mod.IfeChains.ife_id, mod.IfeChains.index)

            if query.count() == 0:
                self.logger.warn("No ifes found for %s" % pdb)
                return []

            grouped = it.imap(result2dict, query)
            grouped = it.groupby(grouped, op.itemgetter('id'))
            groups = []
            for ife_id, chains in grouped:
                chains = list(chains)
                groups.append({
                    'id':  ife_id,
                    'pdb': chains[0]['pdb'],
                    'bp': chains[0]['bp'],
                    'name': chains[0]['name'],
                    'length': chains[0]['length'],
                    'species': chains[0]['species'],
                    'chains': chains,
                    'db_id': chains[0]['db_id'],
                    'resolution': chains[0]['resolution'],
                    'method': chains[0]['method'],
                })

        if not groups:
            raise core.InvalidState("No stored ifes for %s", pdb)

        self.logger.info("Found %i ifes for %s", len(groups), pdb)
        return groups

    def discrepancies(self, groups):
        """Load the discrepancies for the given groups. If use_discrepancy is
        False this will return an empty dictionary. The returned data structure
        will be a dictionary of dictionaries where the final values are Bools.
        The keys in each dictionary are the chain ids which have been aligned.

        :param list groups: The list of groups to use.
        :returns: A nested dictionary of dictionaries.
        """

        if not self.use_discrepancy:
            return {}

        chain_ids = []
        for group in groups:
            chain_ids.append(group['db_id'])

        with self.session() as session:
            sim = mod.ChainChainSimilarity
            query = session.query(sim).\
                filter(sim.chain_id_1.in_(chain_ids)).\
                filter(sim.chain_id_2.in_(chain_ids))

            discrepancy = coll.defaultdict(dict)
            for result in query:
                id1 = result.chain_id_1
                id2 = result.chain_id_2
                discrepancy[id1][id2] = result.discrepancy
                discrepancy[id2][id1] = result.discrepancy

        discrepancy = dict(discrepancy)
        if not discrepancy and self.use_discrepancy:
            raise core.InvalidState("No discrepancy data to cluster with")
        return discrepancy

    def alignments(self, chains):
        """Load the data about alignments between all pairs of chains.

        :chains: A list of chains.
        :returns: A dictionary of dictionaries.
        """

        pdbs = list(set(chain['pdb'] for chain in chains))
        helper = cr.Helper(self.config, self.session.maker)
        alignments = helper.aligned_chains(pdbs, good=True)
        if not alignments:
            raise core.InvalidState("No alignments loaded")
        return alignments

    def has_good_alignment(self, all_alignments, group1, group2):
        """Check if the longest chain of the two groups align well.

        If there no alignments for the first group, or if there is no alignment
        between the first and second we return False. Otherwise we return True
        if there is a good alignment between the two groups.

        :param dict all_alignments: The dictionary of alignments
        :param dict group1: The first group.
        :param dict group2: The second group.
        :returns: True if there is a good alignment.
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
        """Check if the longest chains of the two groups agree. Species are
        similar if they are the same, or if either one is either synthetic or
        None.

        :param dict group1: The first group.
        :param dict group2: The second group.
        :returns: True if the two species are similar.
        """

        if not self.use_species:
            self.logger.debug("Ignoring species assignments")
            return True

        species1 = group1['species']
        species2 = group2['species']
        if species1 != SYNTHEIC[0] and species2 != SYNTHEIC[0] and \
                species1 is not None and species2 is not None \
                and species1 != species2:
            self.logger.debug("Splitting %s, %s: Different species (%i, %i)",
                              group1['id'], group2['id'], species1, species2)
            return False

        self.logger.debug("Good species: %s, %s", group1['id'], group2['id'])
        return True

    def has_good_discrepancy(self, all_discrepancy, group1, group2):
        """Check if two groups have a good discrepancy or not. The discrepancy
        matrix should be a nested dictionary as produced by the `discrepancy`
        method. This only uses the longest chain of each group.

        Because discrepancies were not always computed, there are some special
        cases to consider. If use_discrepancy is False, or if not discrepancies
        are computed then this always returns True.  If either chain is not
        valid for discrepancy as defined by
        `pymotifs.utils.discrepancy.valid_chain` then this will return True.

        If no discrepancies are computed for the first then this returns True.

        Otherwise, if no discrepancy has been computed between two the two
        groups then this returns False. We return true if the discrepancy is
        less than NR_DISCREPANCY_CUTOFF.

        :param dict all_discrepancy: The dictionary of all discrepancies.
        :param dict group1: The first group.
        :param dict group2: The second group.
        :returns: True if there is a good discrepancy between the two groups.
        """

        if not self.use_discrepancy:
            return True

        db_id1 = group1['db_id']
        db_id2 = group2['db_id']
        if not disc.valid_chain(group1) or not disc.valid_chain(group2):
            self.logger.debug("Chains not valid for discrepancy")
            return True

        if not all_discrepancy:
            self.logger.debug("No discrepancies, so ignoring")
            return True

        if db_id1 not in all_discrepancy:
            self.logger.warning("No computed discrepancy for %s", group1['id'])
            return True

        discrepancy = all_discrepancy[db_id1]
        if db_id2 not in discrepancy:
            self.logger.debug("Splitting %s %s: No discrepancy between them",
                              group1['id'], group2['id'])
            return False

        if discrepancy[db_id2] > NR_DISCREPANCY_CUTOFF or discrepancy[db_id2] < 0:
            self.logger.debug("Splitting %s %s: Too large discrepancy %f",
                              group1['id'], group2['id'], discrepancy[db_id2])
            return False

        self.logger.debug("Good discrepancy %s %s", group1['id'], group2['id'])
        return True

    def is_hard_coded_join(self, group1, group2):
        """Check if two groups have a hard coded joining. This uses
        EQUIVALENT_PAIRS from the constants. This should be a set of pairs of
        semi IFE is like ('pdb1|chain1', 'pdb2|chain2'). Note the lack of model
        in the semi IFE id.

        :param dict group1: The first group.
        :param dict group2: The second group.
        :returns: True if there has been a hard coded join between two groups.
        """

        def simple(ife_id):
            parts = ife_id.split('|')
            return (parts[0], parts[-1])

        id1 = (simple(group1['id']), simple(group2['id']))
        id2 = (simple(group2['id']), simple(group1['id']))
        return id1 in EQUIVALENT_PAIRS or id2 in EQUIVALENT_PAIRS

    def are_equivalent(self, alignments, discrepancies, group1, group2):
        """Determine if two chains are equivalent. This requires that the
        chains align well and are not of conflicting species.

        :param dict alignments: A dictionary of dictionaries about the
        alignment between chains.
        :param dict discrepancies: The alignments dictionaries.
        :param dict group1: The first chain.
        :param dict group2: The second chain.
        :returns: True if the two chains are equivalent, False otherwise
        """

        if group1 == group2:
            return True

        if self.is_hard_coded_join(group1, group2):
            return True

        return self.are_similar_species(group1, group2) \
            and self.has_good_alignment(alignments, group1, group2) \
            and self.has_good_discrepancy(discrepancies, group1, group2)

    def validate(self, connections, group):
        """This validates all groups to check if all pairs are connected.

        :param dict connections: All pairs of connections
        :param list group: The group to validate.
        """

        pairs = it.product(group, group)
        pairs = it.ifilter(lambda (a, b): a != b, pairs)
        pairs = it.ifilter(lambda (a, b): b not in connections[a], pairs)
        pairs = list(pairs)
        for pair in pairs:
            self.logger.debug("Pair %s, %s not connected", *pair)
        self.logger.debug("%i pairs are not connected", len(pairs))

    def split_by_species(self, group):
        """Split a group by species. This will put all members with the same
        species into their own group. The ones with synthetic or unknown are
        placed with the largest group.

        :param list group: A list of chains to split by species.
        :returns: The list of species split by group.
        """

        species = coll.defaultdict(list)
        for entry in group:
            name = entry['species']
            if name is None or name == SYNTHENIC_SPECIES_ID:
                name = SYNTHENIC_SPECIES_ID
            species[name].append(entry)

        self.logger.debug("Found groups based on species: %s",
                          ', '.join(str(s) for s in species.keys()))
        species = dict(species)
        unknown = species.pop(SYNTHENIC_SPECIES_ID, [])
        if not species:
            return [group]

        ids = sorted(species.keys(), key=lambda k: (len(species[k]), k))
        largest_group = max(ids, key=lambda k: len(species[k]))
        species[largest_group].extend(unknown)

        merge = lambda group: ', '.join(g['id'] for g in group)
        self.logger.debug("Produced %i groups: %s", len(species),
                          '; '.join(merge(g) for g in species.values()))

        return species.values()

    def enforce_species_splitting(self, group):
        """For groups over the NR_MIN_HOMOGENEOUS_SIZE we require that all
        members of the group have the same species (excluding unknown or
        synthetic species of course). We enforce this requirement here. The
        splitting is done in spilt_by_species. If must_enforce_single_species
        is False then this will simply return the given group.

        :param list group: List of chains in a single group.
        :returns: A list of each subgroup with consistent species, in no
        particular order.
        """

        if not self.must_enforce_single_species:
            return [group]

        max_length = max(group, key=op.itemgetter('length'))
        max_length = max_length['length']
        if max_length < NR_MIN_HOMOGENEOUS_SIZE or len(group) == 1:
            return [group]

        self.logger.debug("Enforcing species splitting of %s",
                          ', '.join(c['id'] for c in group))
        return self.split_by_species(group)

    def missing_ifes(self, matrix, ifes):
        """Check if there are any missing

        :param dict matrix: A dictonary where
        :returns: The list of db_ids whic are not in the matrix.
        """

        missing = it.imap(op.itemgetter('db_id'), ifes)
        missing = it.ifilterfalse(lambda c: c in matrix, missing)
        return bool(list(missing))

    def pairs(self, chains, alignments, discrepancies):
        """Generate an iterator of all equivalent pairs of chains.

        :chains: The chains to build pairs of.
        :alignments: Alignments to use for checking validity.
        :discrepancies: Discrepancies of the chains to use.
        :returns: An iterable of all valid pairs.
        """

        # It should be possible to rewrite this section such that it only
        # iterates over the chains which have been aligned. In addition, it is
        # possible to change the alignment loading to only load the good
        # alignments. If both these changes are done the run time is decreased
        # greatly (I think to ~30 min). For some reason this produces the
        # incorrect results though. This is likely been due to bugs that have
        # been fixed now.
        equiv = ft.partial(self.are_equivalent, alignments, discrepancies)
        pairs = it.combinations(chains, 2)
        return it.ifilter(lambda p: equiv(*p), pairs)

    def connections(self, chains, alignments, discrepancies):
        """Create a graph connections between all chains.

        :param list chains: This chains to build connections with.
        :param dict alignments: The alignments to use.
        :param dict discrepancies: The discrepancies to use.
        :returns: A dictionary of connections to use.
        """

        graph = coll.defaultdict(set)
        for chain in chains:
            graph[chain['id']].add(chain['id'])

        for chain1, chain2 in self.pairs(chains, alignments, discrepancies):
            self.logger.debug("Equivalent: %s %s", chain1['id'], chain2['id'])
            graph[chain1['id']].add(chain2['id'])
            graph[chain2['id']].add(chain1['id'])

        self.logger.info("Created %i connections", len(graph))
        return graph

    def build_groups(self, graph):
        """Group the ifes based upon the given graph. This will group the ifes
        based upon transitivity.

        :param dict graph: A dictionary of id -> set(id...)
        :returns: A list of lists where the inner list is a list of the ife ids
        in a group.
        """

        return cs.find_connected(graph).values()

    def group(self, chains, alignments, discrepancies):
        """Group all chains into connected components.

        :param list chains: A list of chains to group.
        :param dict alignments: The loaded alignments.
        :param dict discrepancies: The loaded discrepancies.
        :returns: A list of lists of the connected components.
        """

        mapping = {}
        for chain in chains:
            mapping[chain['id']] = chain

        groups = []
        graph = self.connections(chains, alignments, discrepancies)
        for ids in self.build_groups(graph):
            self.validate(graph, list(ids))
            group = [mapping[id] for id in ids]
            for subgroup in self.enforce_species_splitting(group):
                groups.append(subgroup)

        return groups

    def all_ifes(self, pdbs):
        ifes = it.imap(self.ifes, pdbs)
        ifes = it.chain.from_iterable(ifes)
        ifes = it.ifilter(self.valid_ife, ifes)
        ifes = list(ifes)
        if not ifes:
            raise core.InvalidState("No ifes found in given pdbs")
        return ifes

    def __call__(self, pdbs, **kwargs):
        """Group all chains in the given list of pdbs to form equivalence
        classes.

        :param list pdbs: A list of pdb ids to group the chains in.
        :returns: A list of groups dictionaries. These will contain one entry,
        'members' a list of chains in this group. These chains will be sorted
        by quality as defined by ranking_key. The chain data structure will
        be a dictionary as returned by merged_chains, but with an extra key
        'rank' which indicates the rank of the chain.
        """

        self.logger.info("Will group %i pdbs", len(pdbs))
        ifes = self.all_ifes(pdbs)
        self.logger.info("Found %i ifes to cluster", len(ifes))
        alignments = self.alignments(ifes)
        discrepancy = self.discrepancies(ifes)

        grouped = set()
        groups = []
        for group in self.group(ifes, alignments, discrepancy):
            members = []
            sorted_group = sorted(group, key=ranking_key)
            for index, member in enumerate(sorted_group):
                grouped.add(member['id'])
                member['rank'] = index
                members.append(member)

            groups.append({'members': members})

        ife_ids = set(ife['id'] for ife in ifes)
        if ife_ids != grouped:
            msg = ("Only {grouped} of {total} ifes were grouped, "
                   "missing: {missing}")
            raise core.InvalidState(msg.format(
                grouped=len(grouped),
                total=len(ife_ids),
                missing=sorted(ife_ids - grouped),
            ))

        self.logger.info("Created %i groups with %i ifes",
                         len(groups), len(grouped))

        return groups
