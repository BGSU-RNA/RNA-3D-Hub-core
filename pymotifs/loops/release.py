import datetime

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import releases as rel

from pymotifs.mat_files import Loader as MatLoader


class Loader(core.MassLoader):
    dependencies = set([MatLoader])

    def current_id(self):
        with self.session() as session:
            query = session.query(mod.LoopReleases.loop_releases_id).\
                order_by(desc(mod.LoopReleases.date)).\
                limit(1)

            if query.count() == 0:
                return '0.0'

            return query.one().loop_releases_id

    def data(self, *args, **kwargs):
        now = datetime.datetime.now()
        current = self.current_id()
        next = rel.next_id(current, mode=self.config['release_mode']['loops'])
        return mod.LoopReleases(loop_releases_id=next, date=now)
