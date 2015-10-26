from pymotifs.core import SimpleLoader
from pymotifs.models import IfeInfo

from pymotifs.ife.grouper import Grouper

from pymotifs.chains.loader import Loader as ChainLoader
from pymotifs.interactions.loader import Loader as InteractionLoader


class Loader(SimpleLoader):
    dependencies = set([ChainLoader, InteractionLoader])
    merge_data = True

    def query(self, session, pdb_id):
        return session.query(IfeInfo).filter_by(pdb_id=pdb_id)

    def as_group(self, group):
        integral_count = None
        accompanying_chain_count = None

        return IfeInfo(
            ife_id=group['id'],
            pdb_id=group['pdb'],
            has_integral=bool(integral_count),
            has_accompanying=bool(accompanying_chain_count),
            has_structured=group['chains'][0]['autonomous'],
            chain_count=len(group['chains']),
            integral_chain_count=None,
            accompanying_chain_count=None,
            structured_chain_count=None,
            length=group['summary']['exp_length'],
            bp_count=group['summary']['bp']
        )

    def data(self, pdb_id):
        grouper = Grouper(self.config, self.session.maker)
        return [self.as_group(g) for g in grouper(pdb_id)]
