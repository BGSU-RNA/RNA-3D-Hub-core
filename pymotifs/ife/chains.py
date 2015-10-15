from pymotifs import core
from pymotifs.models import IfeChains
from pymotifs.models import IfeInfo

from pymotifs.ife.grouper import Grouper

from pymotifs.ife.info import Loader as InfoLoader


class Loader(core.Loader):
    dependencies = set([InfoLoader])

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(IfeChains).\
                join(IfeInfo,
                     IfeInfo.ife_id == IfeChains.ife_id).\
                filter(IfeInfo.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(IfeChains.ife_chain_id).\
                join(IfeInfo,
                     IfeInfo.ife_id == IfeChains.ife_id).\
                filter(IfeInfo.pdb_id == pdb)
            ids = [result.ife_chain_id for result in query]

        if not ids:
            self.logger.warning("Nothing to delete for %s", pdb)
            return None

        with self.session() as session:
            query = session.query(IfeChains).\
                filter(IfeChains.ife_chain_id.in_(ids)).\
                delete(synchronize_session=False)

    def data(self, pdb_id, **kwargs):
        grouper = Grouper(self.config, self.session.maker)
        data = []
        for group in grouper(pdb_id):
            for index, chain in enumerate(group['chains']):
                reference = index == 0
                accompanying = not reference and len(group['chains']) > 1
                data.append(IfeChains(
                    chain_id=chain['db_id'],
                    ife_id=group['id'],
                    is_integral=chain['autonomous'],
                    is_structure=chain['structured'],
                    is_accompanying=accompanying
                ))
        return data
