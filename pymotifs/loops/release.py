import datetime as dt

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import releases as rel

from pymotifs.mat_files import Loader as MatLoader


class Loader(core.MassLoader):
    dependencies = set([MatLoader])
    update_gap = dt.timedelta(7)

    def has_data(self, *args, **kwargs):
        """This will always return True because we only want to update if the time
        difference has been large enough.
        """
        return True

    def remove(self, *args, **kwargs):
        self.logger.info("Will never automatically delete loop releases"
                         " or cached data")

    def current_id(self):
        with self.session() as session:
            query = session.query(mod.LoopReleases.loop_release_id).\
                order_by(desc(mod.LoopReleases.date)).\
                limit(1)

            if query.count() == 0:
                return '0.0'

            return query.one().loop_release_id

    def data(self, *args, **kwargs):
        now = dt.datetime.now()
        current = self.current_id()
        next = rel.next_id(current, mode=self.config['release_mode']['loops'])
        return mod.LoopReleases(loop_release_id=next, date=now)
