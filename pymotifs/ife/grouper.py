"""
This module contains a class which groups chains into integrated functional
elements (IFEs).
These elements are the building blocks of equivalence classes and representative sets.
Some chains make an IFE by themselves and that is OK.
When two chains are internally structured, as with eukaryotic LSU and 5.8S, they are grouped together.
When two chains are unstructured, as with two chains making a simple helix, they are grouped together.
tRNAs are structured but mRNAs generally are not in 3D structures, and it will be better not
to make an IFE for a tRNA-mRNA pair.
Hammerhead ribozymes sometimes have one structured chain and one unstructured.
Special rules are needed there.
"""

import itertools as it
import collections as coll

from pymotifs import core
from pymotifs.constants import IFE_EXTERNAL_INTERNAL_FRACTION as CUTOFF
from pymotifs.ife.helpers import IfeLoader
from pymotifs.ife.helpers import IfeGroup


class Grouper(core.Base):
    """
    This is a class to take a list of chains and find all that are
    structured. The main entry point is to simply call it. Other methods are
    available but these only do part of the tasks required for creating
    structured groupings.
    """

    def partition_interactions(self, ifes, interactions):
        """
        Partition the interactions into those that occur between chains
        which are both structured or both unstructured, versus those that occur
        across chains that have differing structured or unstructured types.
        """

        same = coll.defaultdict(dict)
        rest = coll.defaultdict(dict)
        for ife1, ife2 in it.product(ifes, repeat=2):
            chain1 = ife1.chain
            chain2 = ife2.chain
            count = interactions.get(chain1, {}).get(chain2, 0)
            same[chain1][chain2] = 0
            rest[chain1][chain2] = 0
            if ife1.is_structured == ife2.is_structured:
                same[chain1][chain2] = count
            else:
                rest[chain1][chain2] = count
        return dict(same), dict(rest)


    def should_join(self, ife1, ife2, count):
        """
        Detect if two ifes should be joined. This will consider the types of
        interactions and the ratio of internal base pairs to external base
        pairs.

        :ife1: The first ife chain.
        :ife2: The second ife chain.
        :count: The number of external base pairs between these chains.
        :returns: True if the chains should be joined.
        """

        if ife1 == ife2 or not count:
            return False

        # new rule, there must be more than 3 basepairs between chains to be in an IFE
        if count <= 3:
            return False

        if not ife1.is_structured or not ife2.is_structured:
            return True
        count = float(count)  # avoid integer division


        # avoid division by zero in the next line
        return max(count/(0.001+ife1.internal), count/(0.001+ife2.internal)) >= CUTOFF


    def joinable(self, interactions, *ifes):
        """
        This provides an generator which yields all pairs of chains which
        can be joined according to the given interactions and the properties of
        the chains. This will evaluate all possible pairs of ifes to produce
        the ones which can be joined. The logic for detecting which can be
        joined is in `should_join`.

        This is the same code as was previously used, but only ever
        applied to cases where the type of the chains was the same.

        :inters: A dictionary of dictionaries of counts. For details see
        :*ifes: Lists of ifes to find pairs between.
        """

        kwargs = {}
        if len(ifes) == 1:
            kwargs['repeat'] = 2

        for ife1, ife2 in it.product(*ifes, **kwargs):
            count = interactions[ife1.chain][ife2.chain]
            count2 = interactions[ife2.chain][ife1.chain]
            if count != count2:
                msg = "Reflexive counts not equal: %s (%i), %s (%i)"
                self.logger.warning(msg, ife1, count, ife2, count2)
            count = max(count, count2)

            if self.should_join(ife1, ife2, count):
                yield ife1, ife2


    def joinable_s_u(self, interactions, *ifes):
        """
        This provides a generator which yields all pairs of chains which
        can be joined according to the given interactions and the properties of
        the chains. This will evaluate all possible pairs of ifes to produce
        the ones which can be joined. The logic for detecting which can be
        joined is in `should_join`.

        This logic is used specifically in the case some chains are
        structured and some are not.

        :inters: A dictionary of dictionaries of counts. For details see
        :*ifes: Lists of ifes to find pairs between.
        """

        kwargs = {}
        if len(ifes) == 1:
            kwargs['repeat'] = 2

        longest_length = 0
        for ife1, ife2 in it.product(*ifes, **kwargs):
            if ife1.length > longest_length:
                longest_length = ife1.length

        for ife1, ife2 in it.product(*ifes, **kwargs):
            count = interactions[ife1.chain][ife2.chain]
            count2 = interactions[ife2.chain][ife1.chain]
            if count != count2:
                msg = "Reflexive counts not equal: %s (%i), %s (%i)"
                self.logger.warning(msg, ife1, count, ife2, count2)
            count = max(count, count2)

            # treat structured + unstructured differently since it's a new case
            if ife1 == ife2 or not count:
                pass
            # exclude tRNA and require more interactions than between tRNA and mRNA
            elif longest_length < 60 and count > 3:
                # more than 3 basepairs between structured and unstructured
                yield ife1, ife2


    def group(self, chains, inters):
        """
        Group chains from the same pdb into structured units. This will place
        structured chains into their own groups, but allow a non structured
        chain to be grouped with it if present. Non structured chains that have
        no interacting partners (singleton unstructured chains) are also
        present in the output.

        :ifes: A list of chain objects from nr.chains.Info.load.
        :interactions: A dictionary of interactions of the form produced by
        nr.chains.Info.cross_chain_interactions.
        :returns: A list of ife groups formed from the given ifes and
        interactions.
        """

        # groups is a dictionary that maps chain.id like 1ABC|1|X to IfeGroup object
        groups = {}
        for chain in chains:
            groups[chain.id] = IfeGroup(chain)
            # print("Chain %s bps %s structured %s" % (chain.id,chain.bps,chain.is_structured))
            # print("IFEGroup ordered chains %s" % groups[chain.id].chains())

        # print('chains list:', chains)
        # print('groups dict:', groups)

        # Count interactions between both structured or both unstructured
        # chains (same), and between differently-structured chains (rest)
        same, rest = self.partition_interactions(chains, inters)

        # print('same:', same)
        # print('rest:', rest)

        # If both groups are structured or not structured and the ifes can be
        # joined then they should be merged into one group and both ifes should
        # reflect this.
        # Basically connections between these types of ifes are transitive.
        for chain1, chain2 in self.joinable(same, chains):
            # print('Similarly structured joinable chains %s and %s' % (chain1, chain2))
            current = groups[chain1.id]

            # print('current        before merge: %s' % current)
            # print('current chains before merge: %s' % current.chains())

            # merge creates a group with more chains in it; see helpers.py
            # the group id may or may not show all chains
            current.merge(groups[chain2.id])

            # print('current         after merge: %s' % current)
            # print('current chains  after merge: %s' % current.chains())

            # all of the chains in the two merged groups will now point to current
            linked = set(groups[chain1.id].chains() + groups[chain2.id].chains())
            for chain in linked:
                groups[chain.id] = current

            # print('groups dict after   link: %s' % groups)

        # Here if we are merging an unstructured into a structured we only
        # should update the structured one, the unstructured chain may be part
        # of many different structured ifes.
        # Basically connections are not transitive across types of groups.
        structured = [chain for chain in chains if chain.is_structured]
        unstructured = [chain for chain in chains if not chain.is_structured]

        # self.logger.info('142   structured:', structured)
        # self.logger.info('143 unstructured:', unstructured)

        strip = set()
        for chain1, chain2 in self.joinable_s_u(rest, structured, unstructured):
            # print('Differently structured joinable chains %s and %s' % (chain1, chain2))
            current = groups[chain1.id]

            # print('current        before merge: %s' % current)
            # print('current chains before merge: %s' % current.chains())

            # merge creates a group with more chains in it; see helpers.py
            # however, the group id may or may not show all chains
            current.merge(groups[chain2.id])

            # print('current         after merge: %s' % current)
            # print('current chains  after merge: %s' % current.chains())

            # plan to remove the unstructured chain from the groups
            # that is, do not create an IFE for the unstructured chain
            # now that it is added to an IFE
            strip.add(groups[chain2.id].id)
            # print('strip is now %s' % strip)

        # There will be duplicated groups so we must get only the unique ones.
        # We sort to ensure that the ordering is consistent.
        groups = sorted(set(g for g in groups.values() if g.id not in strip))

        # Some times chains A and B form an IFE, and also A, B, and C form an IFE
        # Exclude the smaller IFE in favor of the larger IFE
        # This is a new rule as of 2024-06-20
        maximal_groups = []
        for group1 in groups:
            chain_set1 = set(c.id for c in group1.chains())
            maximal = True
            for group2 in groups:
                chain_set2 = set(c.id for c in group2.chains())
                if chain_set1 != chain_set2 and chain_set1.issubset(chain_set2):
                    maximal = False
                    self.logger.info('Chain set %s is a strict subset of %s' % (chain_set1, chain_set2))
                    break
            if maximal:
                maximal_groups.append(group1)

        return maximal_groups


    def __call__(self, pdb):
        """
        For the given pdb id, group all chains into structured groups.
        This will group them and ensure that they are consistent.
        It will also merge all chains into one merged chain dictionary.

        :pdb: A list of chain dictionaries to group.
        :returns: A list of grouped chain dictionaries for each structured
        group.
        """

        loader = IfeLoader(self.config, self.session.maker)
        chains, interactions = loader(pdb)
        groups = self.group(chains, interactions)
        if not groups:
            raise core.InvalidState("No ifes found for %s" % pdb)
        self.logger.info("Created %i IFEs for %s" % (len(groups), pdb))
        return groups
