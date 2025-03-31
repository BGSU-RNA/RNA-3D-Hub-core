"""
Load all IFE data for a given set of pdbs into the database.
For each PDB id, group chains together into IFEs.
"""

import itertools as it
from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod

from pymotifs.ife.grouper import Grouper

from pymotifs.units.loader import Loader as UnitInfoLoader
from pymotifs.chains.loader import Loader as ChainLoader
from pymotifs.interactions.loader import Loader as InteractionLoader

from pymotifs.constants import STRUCTURED_BP_COUNT


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader, InteractionLoader, UnitInfoLoader])
    merge_data = True

    def query(self, session, pdb_id):

        return session.query(mod.IfeInfo).\
            filter(mod.IfeInfo.pdb_id == pdb_id).\
            filter(mod.IfeInfo.new_style == True)

    def as_group(self, group):

        # print("group.id               %s" % group.id)
        # print("group.pdb              %s" % group.pdb)
        # print("group.model            %s" % group.model)
        # print("chain_count            %s" % len(group))
        # print("has_structured         %s" % bool(group.is_structured))
        # print("has_integral           %s" % bool(group.integral))
        # print("has_accompanying       %s" % (len(group.chains()) > 1))
        # print("structured_chain_count %s" % len(group.chains(structured=True)))
        # print("length                 %s" % group.length)
        # print("bp_count               %s" % group.bps)

        # some structures have so many chains "stapled" together that they
        # greatly exceed the 100 character limit for the ife_id field
        # we'll reduce the length
        # anyway, the database only allows up to 100 characters, hard to change that
        ife_id = group.group_id()
        while len(ife_id) > 100:
            ife_id = "+".join(ife_id.split("+")[:-1])

        chain_count = len(ife_id.split("+"))

        self.logger.info("ife_id: %s" % ife_id)

        if group.length < 2:
            self.logger.info("Skipping IFE %s because it has fewer than two nucleotides" % ife_id)
            return

        yield mod.IfeInfo(
            ife_id=ife_id,
            pdb_id=group.pdb,
            model=group.model,
            chain_count=chain_count,
            has_structured=bool(group.is_structured),
            has_integral=bool(group.integral),
            has_accompanying=len(group.chains()) > 1,
            structured_chain_count=len(group.chains(structured=True)),  # might exceed chain_count
            length=group.length,
            bp_count=group.bps,
            new_style=True)

        for index, chain in enumerate(group.chains()):
            reference = (index == 0)
            integral = (reference or chain.is_structured)
            accompanying = not integral
            yield mod.IfeChains(
                chain_id=chain.db_id,
                model=chain.model,
                ife_id=ife_id,
                is_integral=integral,
                is_structured=chain.is_structured,
                is_accompanying=accompanying,
                index=index)


    def store_many_chains(self, chain_partners, chain_chain_to_count, pdb_id):
        """
        Used when one chain has many partners, as in DNA origami structures
        """
        # make that chain the IFE and skip the rest of this stage
        chain1 = chain_partners[0][0]
        ife_id = "%s|%s|%s" % (pdb_id, 1, chain1)
        self.logger.info("Chain %s has %d partners, making it the IFE" % (chain1, len(chain_partners[0][1])))
        print("Chain %s has %d partners, making it the IFE" % (chain1, len(chain_partners[0][1])))

        # count total number of units in these chains in the PDB file
        total_length = 0
        with self.session() as session:
            UI = mod.UnitInfo
            query = session.query(UI.unit_id).\
                filter(UI.pdb_id == pdb_id).\
                filter(UI.chain.in_(chain_partners[0][1]))
            total_length = len([row for row in query])

        total_bp_count = 0
        for chain2 in chain_partners[0][1]:
            total_bp_count += chain_chain_to_count[chain1][chain2]

        # map chain_name to chain_id for the ife_chains table
        chain_name_to_id = {}
        with self.session() as session:
            CI = mod.ChainInfo
            query = session.query(CI.chain_id,CI.chain_name).\
                filter(CI.pdb_id == pdb_id)
            for row in query:
                chain_name_to_id[row.chain_name] = row.chain_id

        # record the main chain and its data
        yield mod.IfeInfo(
            ife_id=ife_id,
            pdb_id=pdb_id,
            model=1,
            chain_count=len(chain_partners[0][1]),
            has_structured=(chain_chain_to_count[chain1][chain1] > STRUCTURED_BP_COUNT),
            has_integral=True,
            has_accompanying=True,
            structured_chain_count=len(chain_partners[0][1]),
            length=total_length,
            bp_count=total_bp_count,
            new_style=True)
        # record the main chain and its role
        yield mod.IfeChains(
            chain_id=chain_name_to_id[chain1],
            model=1,
            ife_id=ife_id,
            is_integral=True,
            is_structured=(chain_chain_to_count[chain1][chain1] > STRUCTURED_BP_COUNT),
            is_accompanying=False,
            index=0)
        # record all the other chains as accompanying chains
        counter = 1
        for chain2 in sorted(chain_partners[0][1]):
            if chain2 != chain1:
                yield mod.IfeChains(
                    chain_id=chain_name_to_id[chain2],
                    model=1,
                    ife_id=ife_id,
                    is_integral=False,
                    is_structured=(chain_chain_to_count[chain1][chain1] > STRUCTURED_BP_COUNT),
                    is_accompanying=True,
                    index=counter)
                counter += 1

        print("Only supposed to get here with big origami structures")

        return


    def data(self, pdb_id, **kwargs):

        # query for all basepairs in the structure
        # count cWW basepairs within and between chains
        # join with unit_info to make sure chain_index is not null
        chain_chain_to_count = {}
        with self.session() as session:
            UPI = mod.UnitPairsInteractions2024
            UI1 = aliased(mod.UnitInfo)
            UI2 = aliased(mod.UnitInfo)
            query = session.query(UPI.unit_id_1, UPI.unit_id_2).\
                join(UI1, UPI.unit_id_1 == UI1.unit_id).\
                join(UI2, UPI.unit_id_2 == UI2.unit_id).\
                filter(UPI.pdb_id == pdb_id).\
                filter(UPI.f_lwbp == 'cWW').\
                filter(UPI.program == 'fr3d').\
                filter(UI1.chain_index != None).\
                filter(UI2.chain_index != None)

            for row in query:
                fields1 = row.unit_id_1.split('|')
                # only work with model 1
                if fields1[1] == '1':
                    chain1 = fields1[2]
                    if not chain1 in chain_chain_to_count:
                        chain_chain_to_count[chain1] = {}
                        chain_chain_to_count[chain1][chain1] = 0
                    fields2 = row.unit_id_2.split('|')
                    chain2 = fields2[2]
                    if not chain2 in chain_chain_to_count[chain1]:
                        chain_chain_to_count[chain1][chain2] = 0
                    chain_chain_to_count[chain1][chain2] += 1

        # map each chain to its significant cWW basepairing partners
        chain_to_partners = {}
        for chain1 in chain_chain_to_count:
            for chain2 in chain_chain_to_count[chain1]:
                if chain_chain_to_count[chain1][chain2] > 3:
                    if not chain1 in chain_to_partners:
                        chain_to_partners[chain1] = set([chain1])  # include chain1
                    chain_to_partners[chain1].add(chain2)

        # sort chains and set of partners from most to least partners
        chain_partners = sorted(chain_to_partners.items(), key=lambda x: len(x[1]), reverse=True)

        # make sure the list is not empty ...
        if len(chain_partners) > 0:
            # check for DNA origami structures and if found,
            # completely skip the usual code
            # if the chain with the most partners has more than 10 partners,
            if len(chain_partners[0][1]) > 10:
                # make that chain the IFE and skip the rest of this stage
                groups = self.store_many_chains(chain_partners, chain_chain_to_count, pdb_id)
                return groups

        # traditional code for everything but huge origami structures
        grouper = Grouper(self.config, self.session.maker)
        groups = grouper(pdb_id)
        groups = map(self.as_group, groups)
        groups = it.chain.from_iterable(groups)
        groups = list(groups)

        # print("Number of groups: %d" % len(groups))
        # print(groups)

        return groups
