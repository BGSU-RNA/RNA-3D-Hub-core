"""This module contains a class which groups chains into integrated functional
elements. These elements  are the building blocks of non redudant entires and
thus nr sets.
"""

import itertools as it
import collections as coll

from pymotifs import core
from pymotifs.constants import IFE_EXTERNAL_INTERNAL_FRACTION as CUTOFF
from pymotifs.ife.helpers import IfeLoader
from pymotifs.ife.helpers import IfeGroup


class Grouper(core.Base):
    """This is a class to take a list of chains and find all that are
    structured. The main entry point is to simply call it. Other methods are
    available but these only do part of the tasks required for creating
    structured groupings.
    """

    def should_join(self, ife1, ife2, count):
        """Detect if two ifes should be joined. This will consider the types of
        interactions and the ratio of internal base pairs to external base
        pairs.

        :ife1: The first ife chain.
        :ife2: The second ife chain.
        :count: The number of external base pairs between these chains.
        :returns: True if the chains should be joined.
        """

        if ife1 == ife2 or not count:
            return False
        if not ife1.is_structured or not ife2.is_structured:
            return True
        return max(count/ife1.internal, count/ife2.internal) >= CUTOFF

    def parition_interactions(self, ifes, interactions):
        """Partition the interactions into those that occur between chains
        which are both structured or both unstructured, versus those that occur
        across chains that have differing structured or unstructured types.
        """

        same = coll.defaultdict(dict)
        rest = coll.defaultdict(dict)
        for ife1, ife2 in it.product(ifes, repeat=2):
            chain1 = ife1.chain
            chain2 = ife2.chain
            count = interactions.get(chain1, {}).get(chain2, 0L)
            same[chain1][chain2] = 0L
            rest[chain1][chain2] = 0L
            if ife1.is_structured == ife2.is_structured:
                same[chain1][chain2] = count
            else:
                rest[chain1][chain2] = count
        return dict(same), dict(rest)

    def joinable(self, interactions, *ifes):
        """This provides an generator which yields all pairs of chains which
        can be joined according to the given interactions and the properties of
        the chains. This will evaluate all possible pairs of ifes to produce
        the ones which can be joined. The logic for detecting which can be
        joined is in `should_join`.

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

    def group(self, ifes, inters):
        """Group chains from the same pdb into structured units. This will place
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

        groups = {}
        for ife in ifes:
            groups[ife.id] = IfeGroup(ife)

        same, rest = self.parition_interactions(ifes, inters)

        # If both groups are structured or not structured and the ifes can be
        # joined then they should be merged into one group and both ifes should
        # reflect this. Basically connections between these types of ifes
        # are transitive.
        for ife1, ife2 in self.joinable(same, ifes):
            current = groups[ife1.id]
            current.merge(groups[ife2.id])
            linked = set(groups[ife1.id].chains() + groups[ife2.id].chains())
            for ife in linked:
                groups[ife.id] = current

        # Here if we are merging an unstructured into a structured we only
        # should update the structured one, the unstructured chain may be part
        # of many different structured ifes. Basically connections are not
        # transitive across types of groups.
        structured = [ife for ife in ifes if ife.is_structured]
        unstructured = [ife for ife in ifes if not ife.is_structured]
        strip = set()
        for ife1, ife2 in self.joinable(rest, structured, unstructured):
            groups[ife1.id].merge(groups[ife2.id])
            strip.add(groups[ife2.id].id)

        # There will be duplicated groups so we must get only the unique ones.
        # We sort to ensure that the ordering is consistent.
        return sorted(set(g for g in groups.values() if g.id not in strip))

    def __call__(self, pdb):
        """Group all chains into structured groups. This will group them and
        ensure that they are consistent. It will also merge all chains into one
        merged chain dictionary.

        :pdb: A list of chain dictionaries to group.
        :returns: A list of grouped chain dictionaries for each structured
        group.
        """

        loader = IfeLoader(self.config, self.session.maker)
        chains, interactions = loader(pdb)
        return self.group(chains, interactions)
