"""
Load all IFE data for a given set of pdbs into the database.
For each PDB id, it groups chains together into IFEs.
"""

import itertools as it

from pymotifs import core
from pymotifs import models as mod

from pymotifs.ife.grouper import Grouper

from pymotifs.units.loader import Loader as UnitInfoLoader
from pymotifs.chains.loader import Loader as ChainLoader
from pymotifs.interactions.loader import Loader as InteractionLoader


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

        yield mod.IfeInfo(
            ife_id=group.id,
            pdb_id=group.pdb,
            model=group.model,
            chain_count=len(group),
            has_structured=bool(group.is_structured),
            has_integral=bool(group.integral),
            has_accompanying=len(group.chains()) > 1,
            structured_chain_count=len(group.chains(structured=True)),
            length=group.length,
            bp_count=group.bps,
            new_style=True)

        for index, chain in enumerate(group.chains()):
            print(group.id)
            reference = (index == 0)
            integral = (reference or chain.is_structured)
            accompanying = not integral
            # # print(group.id)
            # # print(333333)
            # with self.session() as session:
            #     query = session.query(mod.IfeChains.ife_id).\
            #             filter(mod.IfeChains.ife_id == group.id) 
            # # print(222222222222)
            # # for row in query:
            # #     print(row.ife_id)
            # # print(not query.count())
            # if not query.count():
            #     # print(1111111111)
            #     # print(group.id)
            yield mod.IfeChains(
                chain_id=chain.db_id,
                model=chain.model,
                ife_id=group.id,
                is_integral=integral,
                is_structured=chain.is_structured,
                is_accompanying=accompanying,
                index=index)

    def data(self, pdb_id, **kwargs):
        grouper = Grouper(self.config, self.session.maker)
        groups = grouper(pdb_id)
        groups = it.imap(self.as_group, groups)
        groups = it.chain.from_iterable(groups)
        groups = list(groups)

        return groups
