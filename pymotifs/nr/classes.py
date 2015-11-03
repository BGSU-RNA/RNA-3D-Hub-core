import datetime as dt

from pymotifs import core
from pymotifs.models import NrClasses

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.utils.releases import Release
from pymotifs.utils import tmp

from pymotifs.chains.loader import Loader as ChainLoader
from pymotifs.interactions.loader import Loader as InteractionLoader
from pymotifs.ife.loader import Loader as IfeLoader


class Loader(core.MassLoader):
    dependencies = set([ReleaseLoader, ChainLoader, InteractionLoader,
                        IfeLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days

    def remove(self, *args, **kwargs):
        self.logger.info("No automatic removal of classes or cached data")

    def has_data(self, pdbs, **kwargs):
        if not tmp.load('nr'):
            raise core.Skip("No cached data")

        helper = Release(self.config, self.session.maker)
        release_id = helper.current('nr')

        with self.session() as session:
            query = session.query(NrClasses).\
                filter_by(nr_release_id=release_id)

            return bool(query.count())

    def data(self, pdbs, **kwargs):
        data = []
        classes = tmp.load('nr')
        for klass in classes:
            data.append(NrClasses(name=klass['name']['full'],
                                  nr_release_id=klass['release'],
                                  handle=klass['name']['handle'],
                                  version=klass['name']['version'],
                                  comment=klass['comment'],
                                  resolution=klass['name']['cutoff']))
        return data
