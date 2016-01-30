import datetime as dt

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod

from pymotifs.utils import releases as rel
from pymotifs.nr.builder import Builder

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.chains.species import Loader as ChainSpeciesLoader
from pymotifs.interactions.loader import Loader as InteractionLoader
from pymotifs.ife.loader import Loader as IfeLoader
from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.chain_chain.loader import Loader as ChainChainLoader


class Loader(core.MassLoader):
    dependencies = set([ChainLoader, ChainSpeciesLoader, InteractionLoader,
                        IfeLoader, CorrespondenceLoader, ChainChainLoader])
    update_gap = dt.timedelta(7)

    def has_data(self, *args, **kwargs):
        """This will always return True because we only want to update if the time
        difference has been large enough.
        """
        return True

    def remove(self, *args, **kwargs):
        self.logger.info("Will never automatically delete nr releases"
                         " or cached data")

    def build(self, pdbs, current_release, next_release, **kwargs):
        builder = Builder(self.config, self.session)
        self.cache('nr', builder(pdbs, current_release, next_release))

    def next_id(self, current):
        return rel.next_id(current, mode=self.config['release_mode']['nrlist'])

    def current_id(self):
        with self.session() as session:
            query = session.query(mod.NrReleases.nr_release_id,
                                  mod.NrReleases.index).\
                order_by(desc(mod.NrReleases.date)).\
                limit(1)

            if query.count() == 0:
                return '0.0', 0

            current = query.one()
            return current.nr_release_id, current.index

    def data(self, pdbs, **kwargs):
        now = kwargs.get('before_date', dt.datetime.now())
        current, index = self.current_id()
        next = self.next_id(current)
        self.build(pdbs, current, next, **kwargs)
        return mod.NrReleases(nr_release_id=next,
                              date=now,
                              parent_nr_release_id=current,
                              description=now.strftime("%Y%m%d"),
                              index=index + 1)
