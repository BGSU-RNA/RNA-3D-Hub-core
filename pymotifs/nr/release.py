import datetime as dt

from pymotifs import core
from pymotifs.models import NrReleases

from pymotifs.utils import tmp
from pymotifs.utils.releases import Release
from pymotifs.nr.builder import Builder

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.interactions.loader import Loader as InteractionLoader


class Loader(core.MassLoader):
    dependencies = set([ChainLoader, InteractionLoader])
    update_gap = dt.timedelta(7)

    def has_data(self, *args, **kwargs):
        """This will always return True because we only want to update if the time
        difference has been large enough.
        """
        return True

    def remove(self, *args, **kwargs):
        tmp.cleanup('nr')
        self.logger.info("Will never automatically delete nr releases")

    def build(self, pdbs, current_release, next_release, **kwargs):
        builder = Builder(self.config, self.session)
        tmp.store('nr', builder(pdbs, current_release, next_release))

    def data(self, pdbs, **kwargs):
        now = dt.datetime.now()
        helper = Release(self.config, self.session)
        current = helper.current('nr')
        next = helper.next(current, mode=self.config['release_mode']['nrlist'])
        self.build(pdbs, current, next, **kwargs)
        return NrReleases(nr_release_id=next, date=now)
