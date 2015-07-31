import datetime as dt

from pymotifs import core
from pymotifs.models import NrReleases
from pymotifs.utils.releases import Release

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.interactions import Loader as InteractionLoader


class Loader(core.MassLoader):
    dependencies = set([ChainLoader, InteractionLoader])
    update_gap = dt.timedelta(7)

    def data(self, *args, **kwargs):
        now = dt.datetime.now()
        helper = Release(self.config, self.session.maker)
        current = helper.current('nr')
        next = helper.next(current, mode=self.config['release_mode']['nrlist'])
        return NrReleases(id=next, date=now)
