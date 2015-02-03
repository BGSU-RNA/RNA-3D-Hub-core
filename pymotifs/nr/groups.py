import itertools as it

import networkx as nx
import networkx.algorithms.components.connected as connected

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils.structures import BasePairQueries as BpHelper
from pymotifs.utils.structures import LONG_RANGE_CUTOFF


class Grouper(object):
    """This is a class to determine the grouping of chains in a pdb. This will
    load all non-redundant chains and then group them by similarity. It
    produces a list of chains which are non-redundant. This does not process
    all chains as it only takes the best of each type from files. However, it
    does make a note about those other chains.
    """

    def __init__(self, config, session_maker):
        self.config = config
        self.session = core.Session(session_maker)
        self.bp_helper = BpHelper(session_maker)

    def best(self, chains):
        """Compute the best chain for a list of chains

        :chains: A list of chains to select the best one of.
        :returns: A single entry from the list that is the best. It adds an
        entry to the chains called 'equivalent', which is all other chains in
        that structure which are equivalent.
        """

        def ordering(chain):
            return (float(chain['bp']) / float(chain['length']),
                    float(chain['length']), chain['name'])

        best = max(chains, ordering)
        best['equivalent'] = list(chains)
        best['equivalent'].remove(best)
        return best

    def rna_chains(self, pdb):
        """Get all RNA containing chains for a given pdb.

        :pdb: The PDB file to get chains for.
        :returns: A list of chain ids for that pdb.
        """
        with self.session() as session:
            query = session.query(mod.ChainInfo.chainId).\
                filter_by(entityMacromoleculeType='Polyribonucleotide (RNA)')
            return [result.chainId for result in query]

    def autonomous(self, chains):
        """This method goes through all chains and detects those that are
        autonomous. When a chain is making no interactions with other chains it
        is kept, if the chain is making important interactions with another
        chain then those will be taken as one autonomous unit and kept.

        :chains: A list of chain dictionaries.
        :returns: A list of merged chain objects.
        """
        pass

    def chains(self, pdb):
        """Load all chains from a given pdb. This will get the RNA chains as
        well as load some interaction data about the chains.

        :pdb: A pdb to get the chains for.
        :returns: A list of dictionaries with data about all chains. The data
        is that which is provided by the info method.
        """
        chains = self.rna_chains(pdb)
        automous_chains = self.autonomous(chains)
        return [self.info(pdb, chain) for chain in automous_chains]

    def entity_mapping(self, pdb):
        data = {}
        with self.session() as session:
            query = session.query(mod.ChainInfo).filter_by(pdb=pdb)
            for result in query:
                data[result.chain] = result.entity_id
        return data

    def nr_chains(self, pdb):
        """Get all chains which are non-redundant from a cif file.
        """
        chains = self.chains(pdb)
        entity_mapping = self.entity_mapping(pdb)
        grouped = it.groupby(chains, lambda c: entity_mapping[c['id']])
        return [self.best(list(chains)) for chain in grouped]

    def alignment_info(self, chains):
        """Load the data about alignments between all pairs of chains.

        :chains: A list of chains.
        :returns: A dictionary of dictionaries.
        """
        pass

    def bps(self, pdb, chain):
        """Determine the total basepairs, internal cWW, external cWW,
        and long range basepairs.

        :pdb: A pdb id to lookup.
        :chain: A chain to use.
        :returns: A dictionary with bp, internal, external and lr entries,
        which stores the respective data.
        """
        helper = self.bp_helper
        return {
            'bp': helper.representative(pdb, chain, count=True),
            'lr': helper.representative(pdb, chain, count=True,
                                        range_cutoff=LONG_RANGE_CUTOFF),
            'internal': helper.representative(pdb, chain, count=True,
                                              family='cWW'),
            'external': helper.cross_chain(pdb, chain, count=True,
                                           family='cWW')
        }

    def info(self, pdb, chain):
        """This loads all information about a chain into a dictionary. This
        will load generic information about a chain, such as resolved, length,
        database id, the source and information about basepairing. The
        basepairing information is loaded from `bps`.

        :pdb: The pdb to search.
        :chain: The chain to search.
        :returns: A dictionary with
        """

        data = {
            'id': '%s|%s' % (pdb, chain),
            'name': chain,
            'pdb': pdb,
        }

        with self.session() as session:
            query = session.query(mod.ChainInfo).\
                filter_by(pdb_id=pdb, chainId=chain)

            result = query.one()
            data['db_id'] = result.id
            data['length'] = result.chainLength
            data['source'] = result.source

        data.update(self.bps(pdb, chain))
        return data

    def is_equivalent(self, alignments, chain1, chain2):
        """Determine if two chains are equivalent.

        :alignments: A dictionary of dictionaries about the alignment between
        chains.
        :chain1: The first chain.
        :chain2: The second chain.
        :returns: True if the two chains are equivalent, False otherwise
        """
        if chain1['source'] != chain2['source']:
            return False
        return True

    def group(self, chains, alignments):
        """Group all chains into connected components.

        :chains: A list of chains to group.
        :returns: A list of lists of the connected components.
        """
        graph = nx.Graph()
        for chain1 in chains:
            for chain2 in chains:
                if self.is_equivalent(alignments, chain1, chain2):
                    graph.add_edge(chain1, chain2)
        return connected.connected_components(graph)

    def representative(self, group):
        """Compute the representative for a group of chains.

        :group: A list of chains.
        :returns: The representative entry from the list.
        """
        pass

    def __call__(self, pdbs):
        chains = self.nr_chains(pdbs)
        alignments = self.alignment_info(chains)
        groups = []
        for group in self.group(chains, alignments):
            groups.append({
                'members': group,
                'representative': self.representative(group)
            })
        return groups
