"""
This module contains a class which groups chains into integrated functional
elements (IFEs).
These elements are the building blocks of equivalence classes and representative sets.
Some chains make an IFE by themselves and that is OK.
When two chains are internally structured, as with eukaryotic LSU and 5.8S, they are grouped together.
When two chains are individually unstructured, as with two chains making a simple helix, they are grouped together.
tRNAs are structured but mRNAs generally are not in 3D structures, and it will be better not
to make an IFE for a tRNA-mRNA pair.
Hammerhead ribozymes sometimes have one structured chain and one unstructured.
Special rules are needed there.
"""

import itertools as it
import collections as coll
from collections import defaultdict
from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs.constants import IFE_EXTERNAL_INTERNAL_FRACTION as CUTOFF
from pymotifs import models as mod

from pymotifs.ife.helpers import IfeLoader
from pymotifs.ife.helpers import IfeGroup
from pymotifs.ife.helpers import IfeChain

from pymotifs.utils import structures as st
from pymotifs import utils as ut


def get_chains_and_counts(self, pdb_id):

    # use a standard procedure to get the nucleic acid chains
    helper = st.Structure(self.session.maker)
    chains = sorted(helper.na_chains(pdb_id))

    # query for all basepairs in the structure
    # find the model with the most basepairs and use that
    # join with unit_info to make sure chain_index is not null
    with self.session() as session:
        UPI = mod.UnitPairsInteractions2024
        UI1 = aliased(mod.UnitInfo)
        UI2 = aliased(mod.UnitInfo)
        query = session.query(UPI.unit_id_1, UPI.unit_id_2, UI1.sym_op, UI2.sym_op).\
            join(UI1, UPI.unit_id_1 == UI1.unit_id).\
            join(UI2, UPI.unit_id_2 == UI2.unit_id).\
            filter(UPI.pdb_id == pdb_id).\
            filter(UPI.f_lwbp == 'cWW').\
            filter(UPI.program == 'fr3d').\
            filter(UI1.chain_index != None).\
            filter(UI2.chain_index != None)

        # find the number of basepairs for each model
        model_to_count = defaultdict(int)
        sym_ops = set()
        for row in query:
            fields = row.unit_id_1.split('|')
            model = fields[1]
            model_to_count[model] += 1
            sym_ops.add(row.sym_op)

        if not model_to_count:
            # if no basepairs, then we have no model information
            with self.session() as session:
                UI = mod.UnitInfo
                mquery = session.query(UI.model).\
                    filter(UI.pdb_id == pdb_id).\
                    distinct()
            for row in mquery:
                model_to_count[row.model] += 1

        if not model_to_count:
            model = '1'
        else:
            model = str(min(model_to_count, key=lambda x: (-model_to_count[x], x)))

        if '1_555' in sym_ops or len(sym_ops) == 0:
            sym_op = '1_555'
        else:
            sym_op = sorted(sym_ops)[0]

        # prepare a place for chain counts
        chain_chain_to_count = {}
        for chain1 in chains:
            chain_chain_to_count[chain1] = {}
            for chain2 in chains:
                chain_chain_to_count[chain1][chain2] = 0

        # count the interactions within and between chains
        for row in query:
            fields1 = row.unit_id_1.split('|')
            if fields1[1] == model:
                # only work with the chosen model
                if len(fields1) < 9 or row.sym_op == sym_op:
                    # skip units with a symmetry operator
                    # we don't join those into IFEs
                    fields2 = row.unit_id_2.split('|')
                    if len(fields2) < 9 or row.sym_op == sym_op:
                        chain1 = fields1[2]
                        chain2 = fields2[2]
                        if chain1 in chains and chain2 in chains:
                            # avoid chains like 5DGF|1|C that is Protein#RNA
                            # but we do include PNA chains
                            # but 7KZL|1|B is all PNA but labeled polypeptide(L), so it's not perfect
                            # also 7UID
                            chain_chain_to_count[chain1][chain2] += 1

        # same chain interactions get counted twice, so fix
        for chain in chains:
            chain_chain_to_count[chain][chain] /= 2

    # get additional information about all of the chains
    with self.session() as session:
        query = session.query(
            mod.ChainInfo.chain_id.label('db_id'),
            mod.ChainInfo.chain_name.label('chain'),
            mod.ChainInfo.pdb_id.label('pdb'),
            mod.ChainInfo.chain_length.label('full_length')).\
            filter_by(pdb_id=pdb_id)

        # loop over chains and create IfeChain objects
        IfeChains = []
        for row in query:
            chain = row.chain
            if chain in chains:

                data = ut.result2dict(row)

                # store the best model here
                data['model'] = model

                # note that full_length is sequence length,
                # not the number of resolved nucleotides
                data['length'] = data['full_length']

                data['internal'] = chain_chain_to_count[chain][chain]
                data['bps'] = chain_chain_to_count[chain][chain]

                IfeChains.append(IfeChain(**data))

    return model, IfeChains, chain_chain_to_count


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

        if ife1 == ife2 or count < 2:
            return False

        # if both chains are short and there are a reasonable number of basepairs, join
        if count >= 2 and ife1.length < 10 and ife2.length < 10:
            return True

        # if there are enough basepairs and at least one chain is unstructured, join
        if count >= 3 and (not ife1.is_structured or not ife2.is_structured):
            return True

        # new rule on 2025-03-06.  Many basepairs, join.  Helps with some LSU-5.8S cases.
        # also helps with 9JM0 chains G and H
        if count >= 6:
            return True

        # when both are unstructured and they don't meet the requirement above, don't join
        if count < 3 and not ife1.is_structured and not ife2.is_structured:
            return False

        count = float(count)  # avoid integer division

        # print('ife1.chain %s ife2.chain %s count %s ife1.internal %s ife2.internal %s' % (ife1.chain,ife2.chain,count,ife1.internal, ife2.internal))

        # this rule was devised with structures smaller than the LSU in mind
        # avoid division by zero
        return max(count/(0.001+ife1.internal), count/(0.001+ife2.internal)) >= CUTOFF


    def joinable(self, interactions, *ifes):
        """
        This provides a generator which yields all pairs of chains which
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

            # print('ife1 %s ife2 %s count %s should_join %s' % (ife1.chain,ife2.chain,count,self.should_join(ife1, ife2, count)))


            # Note: code that makes the ife, when all chains are structured or all are
            # unstructured, will add all chains to the ife, even if part_of_ife is False
            # part_of_ife is only useful when structured and unstructured chains are
            # joined together.
            if self.should_join(ife1, ife2, count):
                if ife1.is_structured:
                    ife1.part_of_ife = True
                if ife2.is_structured:
                    ife2.part_of_ife = True
                if count >= 5:
                    ife1.part_of_ife = True
                    ife2.part_of_ife = True

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

            if ife1 == ife2 or not count:
                pass

            elif ife1.length > 1000:
                # avoid joining long structured chains like SSU with short unstructured
                # chains like mRNA
                # Will mess up some cases, though.
                pass

            # plenty of basepairs between the chains, one is unstructured
            elif count >= 6:
                # print('New joinable_s_u 6 case: %s %s %s' % (ife1, ife2, count))
                ife1.part_of_ife = True
                ife2.part_of_ife = True
                yield ife1, ife2

            # exclude tRNA with mRNA ... why?
            elif count > 3 and longest_length < 60:
                # this should be enough for one chain to be accomanying,
                # so we need to say that they are joinable
                # print('New joinable_s_u 3 case: %s %s %s' % (ife1, ife2, count))
                ife1.part_of_ife = True
                ife2.part_of_ife = True
                yield ife1, ife2

            # all chains here are short and these two have more than 3 cWW interactions
            elif count >= 3:
                # more than 3 basepairs between structured and unstructured
                ife1.part_of_ife = True
                # ife2.part_of_ife is False by default, might have been set to True above
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

        # groups is a dictionary that maps chain.id like 1ABC|1|X to IfeGroup objects
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

        # print('joinable is %s' % sorted(self.joinable(same, chains)))

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

        # print('joinable_s_u is %s' % sorted(self.joinable_s_u(rest, structured, unstructured)))

        strip = set()
        for chain1, chain2 in self.joinable_s_u(rest, structured, unstructured):
            # print('Differently structured joinable chains %s and %s' % (chain1, chain2))
            current = groups[chain1.id]

            # print('current        before merge: %s' % current)
            # print('current id     before merge: %s' % current.group_id())
            # print('current chains before merge: %s' % current.chains())

            # merge creates a group with more chains in it; see helpers.py
            # group_id should not include unstructured chain
            current.merge(groups[chain2.id])

            # print('current         after merge: %s' % current)
            # print('current id      after merge: %s' % current.group_id())
            # print('current chains  after merge: %s' % current.chains())

            # also merge the other direction, for structures like 9CF3 where
            # unstructured joins to unstructured but then never to structured
            # But be careful, because then tRNAs get in the same IFE because
            # they are connected by an mRNA!
            if chain2.part_of_ife:
                groups[chain2.id].merge(groups[chain1.id])

            # maybe this is the right way to identify which chains to avoid
            if not chain2.part_of_ife:
                strip.add(chain2.id)

            # strip.add(groups[chain2.id].id)
            # print('ife.grouper variable strip is now %s' % strip)

        # Avoid making an IFE whose key is an unstructured chain but
        # which is linked to structured chains,
        # because those could link together many structured chains,
        # like when one mRNA binds to many tRNAs
        group_list = []
        for chain_id, group in groups.items():
            if not chain_id in strip:
                group_list.append(group)
                # print('adding group  %s' % group.group_id())
                # print('group chains %s' % group.chains())

        # There will be duplicated groups so we must get only the unique ones.
        # We sort to ensure that the ordering is consistent.
        group_set = set(group_list)

        groups = sorted(group_set, key = lambda x: x.group_id())

        # Sometimes chains A and B form an IFE, and also A, B, and C form an IFE
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

                print('maximal group  %s' % group1.group_id())
                for chain in group1.chains():
                    print('  chain %s .structured %s .part_of_ife %s' % (chain.chain, chain.is_structured, chain.part_of_ife))

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

        # using get_chains_and_counts reduces time from 38 seconds to 2 seconds on 7K00
        self.logger.info("Starting get_chains_and_counts for %s" % pdb)
        model, chains, interactions = get_chains_and_counts(self, pdb)
        # print('Found chains %s' % chains)
        # print('Found interactions %s' % interactions)
        self.logger.info("Found chains %s" % chains)
        self.logger.info("Found interactions %s" % interactions)

        # old loader below, before 3/7/2025
        # it's slow
        # self.logger.info("Starting IfeLoader for %s" % pdb)
        # loader = IfeLoader(self.config, self.session.maker)
        # chains, interactions2 = loader(pdb)

        # print('Found chains %s' % chains)
        # print('Found interactions %s' % interactions2)
        # self.logger.info("Found chains %s" % chains)
        # self.logger.info("Found interactions %s" % interactions2)

        self.logger.info("Grouping %d chains for %s" % (len(chains),pdb))

        groups = self.group(chains, interactions)
        if not groups:
            raise core.InvalidState("No ifes found for %s" % pdb)
        self.logger.info("Created %i IFEs for %s" % (len(groups), pdb))
        return groups
