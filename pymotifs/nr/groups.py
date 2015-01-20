import itertools as it

import networkx as nx
import networkx.algorithms.components.connected as connected

from pymotifs import core
from pymotifs.models import ChainInfo


class Grouper(core.Stage):
    """This is a class to determine the grouping of chains in a pdb. This will
    load all non-redundant chains and then group them by similarity. It
    produces a list of chains which are non-redundant. This does not process
    all chains as it only takes the best of each type from files. However, it
    does make a note about those other chains.
    """

    def __init__(self, config, session_maker):
        self.config = config
        self.session = core.Session(session_maker)

    def best(self, chains):
        """Compute the best chain for a list of chains

        :chains: A list of chains to select the best one of.
        :returns: A single entry from the list that is the best. It adds an
        entry to the chains called 'equivalent', which is all other chains in
        that structure which are equivalent.
        """
        pass

    def rna_chains(self, pdb):
        """Get all RNA containing chains for a given pdb.

        :pdb: The PDB file to get chains for.
        :returns: A list of chain ids for that pdb.
        """
        with self.session() as session:
            query = session.query(ChainInfo.chainId).\
                filter_by(entityMacromoleculeType='Polyribonucleotide (RNA)')
            return [result.chainId for result in query]

    def chains(self, pdb):
        """Load all chains from a given pdb. This will get the RNA chains as
        well as load some interaction data about the chains.

        :pdb: A pdb to get the chains for.
        :returns: A list of dictionaries with data about all chains. The data
        is that which is provided by the info method.
        """
        chains = self.rna_chains(pdb)
        return [self.info(pdb, chain) for chain in chains]

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
        pass

    def info(self, pdb, chain):
        data = {
            'id': '%s|%s' % (pdb, chain),
            'name': chain,
            'pdb': pdb,
            'source': self.source(pdb, chain),
        }

        data.update(self.bps(pdb, chain))
        return data

    def is_equivalent(self, chain1, chain2):
        """Determine if two chains are equivalent.

        :chain1: The first chain.
        :chain2: The second chain.
        :returns: True if the two chains are equivalent, False otherwise
        """
        pass

    def group(self, chains):
        """Group all chains into connected components.

        :chains: A list of chains to group.
        :returns: A list of lists of the connected components.
        """
        graph = nx.Graph()
        for chain1 in chains:
            for chain2 in chains:
                if self.is_equivalent(chain1, chain2):
                    graph.add_edge(chain1, chain2)
        return connected.connected_components(graph)

    def __call__(self, pdbs):
        chains = self.nr_chains(pdbs)
        return self.group(chains)
