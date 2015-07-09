"""This module contains a class which groups chains into autonomous units.
These units are the building blocks of non redudant entires and thus nr sets.
"""

import collections as coll

from pymotifs import core
from pymotifs.utils import connectedsets as cs
from pymotifs.autonomous.info_loader import Info


class Grouper(core.Base):
    """This is a class to take a list of chains and find all that are
    autonomous. The main entry point is to simply call it. Other methods are
    available but these only do part of the tasks required for creating
    autonomous groupings.
    """

    cutoffs = {
        'internal': 0.5,
        'internal_count': 5
    }

    def is_autonomous(self, chain):
        """Determine if a chain is autonomous.

        :chain: A chain dictionary.
        :returns: True or False.
        """

        internal = float(chain['internal'])
        if not internal or internal < self.cutoffs['internal_count']:
            return False
        total = internal + float(chain['external'])
        fraction = internal / total
        return fraction >= self.cutoffs['internal']

    def split_group(self, group, interactions):
        """This splits a group that contains more than one autonomous chain. In
        this case what we want is to separate out the chains so that they
        contain only autonomous chain per group. What we do is to redo the
        grouping, but without the autonomous chains to force them apart. Note
        the chains must already have an autonomous entry.

        :group: A list of chains to separate .
        :interactions: The dictonary of interactions as for self.group.
        :returns: A list of groups from this one group.
        """

        autonomous = [c for c in group if c['autonomous']]
        non_auto = [c for c in group if not c['autonomous']]

        result = []
        for current in autonomous:
            tmp_chains = list(non_auto)
            tmp_chains.append(current)
            result.extend(self.group(tmp_chains, interactions))
        return result

    def group(self, chains, interactions):
        """Group chains from the same pdb into autonomous units. This will place
        autonomous chains into their own groups, but allow a non autonomous
        chain to be grouped with it if present. Non autonomous chains that have
        no interacting partners (singleton unstructured chains) are also
        present in the output.

        :chains: A list of chain objects from nr.chains.Info.load.
        :interactions: A dictionary of interactions of the form produced by
        nr.chains.Info.cross_chain_interactions.
        :returns: A list of lists of the chains in a group. Each list is
        suitable for merging and is sorted by the chain names. In addition, all
        entries are sorted by the first chain name in each group. This makes
        sure the results are stable across runs.
        """

        connections = coll.defaultdict(set)
        mapping = {}
        for chain in chains:
            mapping[chain['id']] = chain
            chain['autonomous'] = self.is_autonomous(chain)
            if chain['autonomous']:
                connections[chain['id']].add(chain['id'])
                self.logger.debug("Chain %s is autonomous", chain['id'])
                continue

            connections[chain['id']].add(chain['id'])
            for chain2 in chains:
                if chain == chain2:
                    continue

                if interactions.get(chain['name'], {}).get(chain2['name']):
                    self.logger.debug("Chain %s has interactions with %s",
                                      chain['id'], chain2['id'])
                    connections[chain['id']].add(chain2['id'])

        grouped = []
        groups = cs.find_connected(dict(connections))
        for group in groups.values():
            group_id = ','.join(group)
            self.logger.debug("Generated group %s from interactions",
                              group_id)
            chains = [mapping[name] for name in group]
            group_autonomous = [c for c in chains if c['autonomous']]

            if len(group_autonomous) > 1:
                self.logger.debug("Partitioning %s due to >1 autonomous chain",
                                  group_id)
                grouped.extend(self.split_group(chains, interactions))
            else:
                grouped.append(sorted(chains, key=lambda c: c['name']))

        return sorted(grouped, key=lambda entry: entry[0]['name'])

    def __call__(self, pdb):
        """Group all chains into autonomous groups. This will group them and
        ensure that they are consistent. It will also merge all chains into one
        merged chain dictionary.

        :pdb: A list of chain dictionaries to group.
        :returns: A list of grouped chain dictionaries for each autonomous
        group.
        """

        info = Info(self.config, self.session.maker)
        chains = info.load_all(pdb)
        interactions = info.cross_chain_interactions(chains)
        groups = self.group(chains, interactions)
        return [info.merge(g) for g in groups]
