"""Store a new NR release and cache and NR grouping.

This will compute the next NR release id and store an entry in the NR release
table about it. It will also use the given PDBs to create a new NR grouping.
This grouping will then be cached for usage in future stages.
"""

import datetime as dt

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.utils import releases as rel
from pymotifs.nr.builder import Builder

from pymotifs.chains.info import Loader as ChainLoader
# from pymotifs.chains.species import Loader as ChainSpeciesLoader
from pymotifs.interactions.loader import Loader as InteractionLoader
from pymotifs.ife.loader import Loader as IfeLoader
from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.chain_chain.loader import Loader as ChainChainLoader
from pymotifs.quality.loader import Loader as QualityLoader
from pymotifs.units.loader import Loader as UnitsLoader


class Loader(core.MassLoader):
    dependencies = set([ChainLoader, InteractionLoader,
                        IfeLoader, CorrespondenceLoader, ChainChainLoader,
                        QualityLoader, UnitsLoader])
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
        self.cache(NR_CACHE_NAME, builder(pdbs, current_release, next_release))

    def next_id(self, current):
        return rel.next_id(current, mode=self.config['release_mode']['nrlist'])

    def current_id(self):
        with self.session() as session:
            query = session.query(mod.NrReleases.nr_release_id,
                                  mod.NrReleases.index).\
                order_by(desc(mod.NrReleases.index)).\
                limit(1)

            if query.count() == 0:
                return '0.0', 0

            current = query.one()
            return current.nr_release_id, current.index

    def data(self, pdbs, **kwargs):

        now = kwargs.get('before_date', dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S'))

        if ":" in now:
            nowstring = dt.datetime.strftime(dt.datetime.strptime(now, "%Y-%m-%d %H:%M:%S").date(), "%Y%m%d")
        else:
#            nowstring = now.strftime("%Y%m%d"),
            nowstring = now.replace("-","")

        current, index = self.current_id()
        next = self.next_id(current)
        parent = current
        if current == '0.0':
            parent = next

        # build the current release
        self.build(pdbs, parent, next, **kwargs)

        # store the release id information
        return mod.NrReleases(nr_release_id=next,
                              date=now,
                              parent_nr_release_id=parent,
                              description = nowstring,
                              index=index + 1)
