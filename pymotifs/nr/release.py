import datetime

from pymotifs import core
from pymotifs.models import NrReleases
from pymotifs.utils.releases import Release


class Loader(core.MassLoader):
    def data(self, *args, **kwargs):
        now = datetime.datetime.now()
        helper = Release(self.config, self.session.maker)
        current = helper.current('nr')
        next = helper.next(current, mode=self.config['release_mode']['nrlist'])
        return NrReleases(id=next, date=now)
