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
        integral_count = None
        accompanying_count = len(group['chains']) - 1
        structured_count = sum(g['autonomous'] for g in group['chains'])

        yield mod.IfeInfo(
            ife_id=group['id'],
            pdb_id=group['pdb'],
            has_integral=bool(integral_count),
            has_accompanying=bool(accompanying_count),
            has_structured=group['chains'][0]['autonomous'],
            chain_count=len(group['chains']),
            integral_chain_count=integral_count,
            accompanying_chain_count=accompanying_count,
            structured_chain_count=structured_count,
            length=group['summary']['exp_length'],
            bp_count=group['summary']['bp'])

        for index, chain in enumerate(group['chains']):
            reference = index == 0
            accompanying = not reference and len(group['chains']) > 1
            yield mod.IfeChains(
                chain_id=chain['db_id'],
                ife_id=group['id'],
                is_integral=chain['autonomous'],
                is_structured=chain.get('structured', False),
                is_accompanying=accompanying)

    def data(self, pdb_id, **kwargs):
        grouper = Grouper(self.config, self.session.maker)
        groups = grouper(pdb_id)
        groups = it.imap(self.as_group, groups)
        groups = it.chain.from_iterable(groups)
        return list(groups)
