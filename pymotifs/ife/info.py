from pymotifs.core import SimpleLoader
from pymotifs.models import IfeInfo

from pymotifs.ife.grouper import Grouper

from pymotifs.chains import Loader as ChainLoader
from pymotifs.interactions import Loader as InteractionLoader


class Loader(SimpleLoader):
    dependencies = set([ChainLoader, InteractionLoader])
    merge_data = True

    def query(self, session, pdb_id):
        return session.query(IfeInfo).filter_by(pdb_id=pdb_id)

    def as_group(self, group):
        return IfeInfo(
            id=group['id'],
            pdb_id=group['pdb'],
            has_structured=group['chains'][0]['autonomous'],
            chain_count=len(group['chains']),
            integral_count=None,
            accompanying_count=None,
            length=group['summary']['exp_length'],
            bps=group['summary']['bp']
        )

    def data(self, pdb_id):
        grouper = Grouper(self.config, self.session.maker)
        return [self.as_group(g) for g in grouper(pdb_id)]
