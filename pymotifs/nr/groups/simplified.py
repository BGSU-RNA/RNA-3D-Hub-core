import itertools as it
import collections as coll

import networkx as nx

from pymotifs import core
from pymotifs.utils import result2dict

from pymotifs.models import ChainInfo
from pymotifs.models import ChainSpecies
from pymotifs.models import AutonomousInfo
from pymotifs.models import AutonomousChains
from pymotifs.models import ChainChainSimilarity

from pymotifs.utils import connectedsets as cs
from pymotifs.utils import correspondence as cr
from pymotifs.utils.structures import SYNTHEIC


def ordering(group):
    rep = group['representative']
    return (rep['summary']['exp_length'], rep['chains'][0]['sequence'])


class Grouper(core.Base):
    """This is a class to determine the grouping of chains in a pdb. This will
    load all non-redundant chains and then group them by similarity. It
    produces a list of chains which are non-redundant. This does not process
    all chains as it only takes the best of each type from files. However, it
    does make a note about those other chains.
    """

    cutoffs = {
        'discrepancy': 0.5
    }

    def chains(self, pdb):
        """Load all autonomous chains from a given pdb. This will get the RNA chains as
        well as load some interaction data about the chains.

        :pdb: A pdb to get the chains for.
        :returns: A list of dictionaries with data about all chains. The data
        is that which is provided by the info method.
        """

        with self.session() as session:
            query = session.query(AutonomousChains.chain_id.label('db_id'),
                                  AutonomousChains.is_reference,
                                  AutonomousChains.is_autonomous,
                                  AutonomousChains.autonomous_id.label('id'),
                                  ChainInfo.pdb_id.label('pdb'),
                                  ChainInfo.chain_length.label('exp_length'),
                                  ChainInfo.sequence,
                                  ChainSpecies.species_id.label('source')).\
                join(AutonomousInfo,
                     AutonomousInfo.id == AutonomousChains.autonomous_id).\
                join(ChainInfo,
                     ChainInfo.id == AutonomousChains.chain_id).\
                join(ChainSpecies,
                     ChainSpecies.chain_id == ChainInfo.id).\
                filter(AutonomousInfo.pdb_id == pdb).\
                order_by(AutonomousChains.autonomous_id)

            grouped = it.groupby(it.imap(result2dict, query),
                                 lambda g: g['id'])
            groups = []
            for group_id, chains in grouped:
                chains = sorted(chains, key=lambda c: c['is_reference'])
                groups.append({
                    'id':  group_id,
                    'pdb': chains[0]['pdb'],
                    'chains': chains
                })

        return groups

    def discrepancy(self, groups):
        chain_ids = []
        for group in groups:
            chain_ids.extend(chain['db_id'] for chain in group['chains'])

        with self.session() as session:
            query = session.query(ChainChainSimilarity).\
                filter(ChainChainSimilarity.chain_id1.in_(chain_ids)).\
                filter(ChainChainSimilarity.chain_id2.in_(chain_ids))

            discrepancy = coll.defaultdict(dict)
            for result in query:
                id1 = result.chain_id1
                id2 = result.chain_id2
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

        db_id1 = group1['chains'][0]['db_id']
        db_id2 = group2['chains'][0]['db_id']
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

        source1 = group1['chains'][0]['source']
        source2 = group2['chains'][0]['source']
        if source1 != SYNTHEIC[0] and source2 != SYNTHEIC[0] and \
                source1 is not None and source2 is not None \
                and source1 != source2:
            self.logger.debug("Splitting %s, %s: Different sources (%i, %i)",
                              group1['id'], group2['id'], source1, source2)
            return False

        self.logger.debug("Good species: %s, %s", group1['id'], group2['id'])
        return True

    def has_good_discrepancy(self, all_discrepancy, group1, group2):
        db_id1 = group1['chains'][0]['db_id']
        db_id2 = group2['chains'][0]['db_id']
        if db_id1 not in all_discrepancy:
            self.logger.debug("Splitting %s %s: No discrepancy for %s",
                              group1['id'], group2['id'], group1['id'])
            return False

        discrepancy = all_discrepancy[db_id1]
        if db_id2 not in discrepancy:
            self.logger.debug("Splitting %s %s: No discrepancy between them",
                              group1['id'], group2['id'])
            return False

        if discrepancy[db_id2] > self.cutoffs['discrepancy']:
            self.logger.debug("Splitting %s %s: Too large discrepancy %f",
                              group1['id'], group2['id'], discrepancy[db_id2])
            return False

        self.logger.debug("Good discrepancy: %s, %s",
                          group1['id'], group2['id'])
        return True

    def are_equivalent(self, alignments, discrepancies, group1, group2):
        """Determine if two chains are equivalent. This requires that the
        chains align well and are not of conflicting species.

        :alignments: A dictionary of dictionaries about the alignment between
        chains.
        :chain1: The first chain.
        :chain2: The second chain.
        :returns: True if the two chains are equivalent, False otherwise
        """

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

    def group(self, chains, alignments, discrepancies):
        """Group all chains into connected components.

        :chains: A list of chains to group.
        :returns: A list of lists of the connected components.
        """

        graph = coll.defaultdict(set)
        G = nx.Graph()
        G.add_nodes_from(c['id'] for c in chains)
        mapping = {}
        for chain1 in chains:
            mapping[chain1['id']] = chain1

            for chain2 in chains:
                if chain1 == chain2 or \
                        self.are_equivalent(alignments, discrepancies, chain1,
                                            chain2):
                    self.logger.debug("Equivalent: %s %s",
                                      chain1['id'], chain2['id'])

                    graph[chain1['id']].add(chain2['id'])
                    graph[chain2['id']].add(chain1['id'])
                    G.add_edge(chain1['id'], chain2['id'])
                else:
                    self.logger.debug("Non-equivalent: %s %s",
                                      chain1['id'], chain2['id'])

        # shortest = nx.all_shortest_paths(G, source='2QNH|2', target='2B64|V')
        # print([p for p in shortest])
        connected = cs.find_connected(graph).values()
        # conneced = nx.algorithms.clique.find_cliques(G)

        groups = []
        for ids in connected:
            self.validate(graph, list(ids))
            groups.append(mapping[id] for id in ids)

        return groups

    def representative(self, chains):
        """Compute the representative for a group of chains.

        :group: A list of chains.
        :returns: The representative entry from the list.
        """

        return chains[0]

    def __call__(self, pdbs, **kwargs):
        """Group all chains in the given list of pdbs.

        :pdbs: A list of pdb ids to group the chains in.
        :returns: A list of groups dictonaries. These will contain two entries,
        'representative' the representative chain of that group and 'members' a
        list of chains in this group. The chain data structure will be a
        dictonary as returned by merged_chains.
        """

        chains = it.imap(self.chains, pdbs)
        chains = list(it.chain.from_iterable(chains))

        groups = []
        alignments = self.alignments(chains)
        discrepancy = self.discrepancy(chains)

        for group in self.group(chains, alignments, discrepancy):
            members = list(group)
            representative = self.representative(members)
            members.remove(representative)
            groups.append({
                'members': members,
                'representative': representative
            })

        return groups
