from pymotifs.core import SimpleLoader
from pymotifs.models import AutonomousGroups

from pymotifs.autonomous.grouper import Grouper


class Loader(SimpleLoader):

    def query(self, session, pdb_id):
        return session.query(AutonomousGroups).filter_by(pdb_id=pdb_id)

    def as_group(self, group):
        return AutonomousGroup(
            id=group['id'],
            pdb_id=group['pdb'],
            has_autonomous=group['chains'][0]['autonomous'],
            chain_count=len(group['chains']),
            length=group['exp_length']
        )

    def data(self, pdb_id):
        grouper = Grouper(self.config, self.session.maker)
        return [self.as_group(g) for g in grouper(pdb_id)]
