from pymotifs.core import SimpleLoader
from pymotifs.models import AutonomousInfo

from pymotifs.autonomous.grouper import Grouper

from pymotifs.chains import Loader as ChainLoader
from pymotifs.interactions import Loader as InteractionLoader


class Loader(SimpleLoader):
    dependencies = set([ChainLoader, InteractionLoader])

    def query(self, session, pdb_id):
        return session.query(AutonomousInfo).filter_by(pdb_id=pdb_id)

    def as_group(self, group):
        return AutonomousInfo(
            id=group['id'],
            pdb_id=group['pdb'],
            has_autonomous=group['chains'][0]['autonomous'],
            chain_count=len(group['chains']),
            length=group['summary']['exp_length']
        )

    def data(self, pdb_id):
        grouper = Grouper(self.config, self.session.maker)
        return [self.as_group(g) for g in grouper(pdb_id)]
