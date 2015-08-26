from pymotifs import core
from pymotifs.models import AutonomousChains
from pymotifs.models import AutonomousInfo

from pymotifs.ife.grouper import Grouper

from pymotifs.ife.info import Loader as InfoLoader


class Loader(core.Loader):
    dependencies = set([InfoLoader])

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(AutonomousChains).\
                join(AutonomousInfo,
                     AutonomousInfo.id == AutonomousChains.autonomous_id).\
                filter(AutonomousInfo.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(AutonomousChains.id).\
                join(AutonomousInfo,
                     AutonomousInfo.id == AutonomousChains.autonomous_id).\
                filter(AutonomousInfo.pdb_id == pdb)
            ids = [result.id for result in query]

        if not ids:
            self.logger.warning("Nothing to delete for %s", pdb)
            return None

        with self.session() as session:
            query = session.query(AutonomousChains).\
                filter(AutonomousChains.id.in_(ids)).\
                delete(synchronize_session=False)

    def data(self, pdb_id, **kwargs):
        grouper = Grouper(self.config, self.session.maker)
        data = []
        for group in grouper(pdb_id):
            for index, chain in enumerate(group['chains']):
                reference = index == 0
                accompanying = not reference and len(group['chains']) > 1
                data.append(AutonomousChains(
                    chain_id=chain['db_id'],
                    autonomous_id=group['id'],
                    is_autonomous=chain['autonomous'],
                    is_reference=reference,
                    is_accompanying=accompanying
                ))
        return data
