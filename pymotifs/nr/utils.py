from pymotifs import core

from pymotifs.constants import NR_CACHE_NAME


class BaseLoader(core.Loader):
    def query(self, session, release):
        if not self.table:
            raise core.InvalidState("Must define table for BaseLoaders")

        return session.query(self.table).\
            filter(self.table.nr_release_id == release)

    def to_process(self, args, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("No precomputed grouping to store")
        return [data['release']]

    def has_data(self, release, **kwargs):
        with self.session() as session:
            query = self.query(session, release)
            return bool(query.count())

    def remove(self, release, **kwargs):
        with self.session() as session:
            query = self.query(session, release)
            query.delete(synchronize_session=False)
        self.logger.warning("No automatic removal of cached NR data")
