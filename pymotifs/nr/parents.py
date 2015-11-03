import datetime as dt

from pymotifs import core

from pymotifs.models import NrClassParents

from pymotifs.utils import tmp
from pymotifs.nr.classes import Loader as ClassLoader


class Loader(core.MassLoader):
    dependencies = set([ClassLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days
    allow_no_data = True

    def has_data(self, *args, **kwargs):
        grouping = tmp.load('nr')
        if not grouping:
            raise core.Skip("No precomputed grouping to store")

        release_id = grouping[0]['release']
        with self.session() as session:
            query = session.query(NrClassParents).\
                filter_by(nr_release_id=release_id)
            return bool(query.count())

    def remove(self, *args, **kwargs):
        self.logger.info("No automatic removal of parents or cached data")

    def parents(self, grouping):
        if not grouping:
            raise core.InvalidState("Must give grouping")

        data = []
        for group in grouping:
            for parent in group['parents']:
                data.append({
                    'nr_class_id': group['class_id'],
                    'release_id': group['release'],
                    'nr_class_parent_id': parent['name']['class_id']
                })

        return data

    def data(self, pdbs, **kwargs):
        grouping = tmp.load('nr')
        for parent in self.parents(grouping):
            yield NrClassParents(**parent)
        tmp.cleanup('nr')
