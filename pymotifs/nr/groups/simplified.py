import itertools as it
import collections as coll
from pprint import pprint

import networkx as nx

from pymotifs import core
from pymotifs.nr import connectedsets as cs
from pymotifs.nr.groups.autonomous import AutonomousGrouper
from pymotifs.utils import correspondence as cr
from pymotifs.nr.chains import Info as ChainLoader
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

    def chains(self, pdb):
        """Load all autonomous chains from a given pdb. This will get the RNA chains as
        well as load some interaction data about the chains.

        :pdb: A pdb to get the chains for.
        :returns: A list of dictionaries with data about all chains. The data
        is that which is provided by the info method.
        """

        loader = ChainLoader(self.config, self.session.maker)
        grouper = AutonomousGrouper(self.config, self.session.maker)
        return grouper(loader.load_all(pdb))[0]

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

    def has_good_discrepancy(self, group1, group2):
        return True

    def are_equivalent(self, alignments, group1, group2):
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
            self.has_good_discrepancy(group1, group2)

    def validate(self, connections, group):
        pairs = it.product(group, group)
        pairs = it.ifilter(lambda (a, b): a != b, pairs)
        pairs = it.ifilter(lambda (a, b): b not in connections[a], pairs)
        pairs = list(pairs)
        for pair in pairs:
            self.logger.debug("Pair %s, %s not connected", *pair)
        self.logger.debug("%i pairs are not connected", len(pairs))

    def group(self, chains, alignments):
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
                        self.are_equivalent(alignments, chain1, chain2):
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

        # def ordering(chain):
        #     if chain['bp'] == 0:
        #         bp_nt = 0
        #     else:
        #         bp_nt = float(chain['bp']) / float(chain['length'])

        #     name = chain['name']
        #     if isinstance(name, list):
        #         name = min(name)

        #     return (bp_nt, chain['length'])

        # ordered = sorted(chains, key=ordering)
        # _, best = next(it.groupby(ordered, key=ordering))
        # best = list(best)

        # # If there is more than 1 chain that is the best, we pick the one
        # # with the smallest name, ie, A instead of C.
        # if len(best) > 1:
        #     return min(best, key=lambda b: b['name'])
        # return best[0]

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

        for group in self.group(chains, alignments):
            members = list(group)
            representative = self.representative(members)
            members.remove(representative)
            groups.append({
                'members': members,
                'representative': representative
            })

        return sorted(groups, key=ordering)
