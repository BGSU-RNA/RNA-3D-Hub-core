"""Load IFE data. This will load all IFE data for a given set of pdbs into the
database.
"""

import itertools as it

from pymotifs import core
from pymotifs import models as mod

from pymotifs.ife.grouper import Grouper

from pymotifs.chains.loader import Loader as ChainLoader
from pymotifs.interactions.loader import Loader as InteractionLoader


class Loader(core.SimpleLoader):
    dependencies = set([ChainLoader, InteractionLoader])
    merge_data = True

    def query(self, session, pdb_id):
        return session.query(mod.IfeInfo).filter_by(pdb_id=pdb_id)

    def as_group(self, group):
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
            reference = (index == 0)
            integral = (reference or chain.is_structured)
            accompanying = not integral
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
        return list(groups)
