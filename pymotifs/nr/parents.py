import datetime as dt

from pymotifs import core

from pymotifs import models as mod

from pymotifs.nr.builder import Builder
from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.classes import Loader as ClassLoader


class Loader(core.MassLoader):
    dependencies = set([ReleaseLoader, ChainLoader, ClassLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days
    allow_no_data = True
    table = mod.NrClassParents

    def has_data(self, *args, **kwargs):
        grouping = self.cached('nr')
        if not grouping:
            raise core.Skip("No precomputed grouping to store")

        release_id = grouping[0]['release']
        with self.session() as session:
            query = session.query(mod.NrClassParents).\
                filter_by(nr_release_id=release_id)
            return bool(query.count())

    def remove(self, *args, **kwargs):
        self.logger.info("No automatic removal of parents or cached data")

    def mapping(self, grouping):
        helper = Builder(self.config, self.session)
        release_id = grouping[0]['release']
        classes = [g['name']['full'] for g in grouping]
        return helper.class_id_mapping(classes, release_id)

    def parents(self, grouping, mapping):
        if not grouping:
            raise core.InvalidState("Must give grouping")

        if not mapping:
            raise core.InvalidState("Must give mapping")

        data = []
        for group in grouping:
            if group['name']['full'] not in mapping:
                raise core.InvalidState("Group %s not in mapping" % group)
            nr_class_id = mapping[group['name']['full']]

            for parent in group['parents']:
                data.append({
                    'nr_class_id': nr_class_id,
                    'nr_release_id': group['release'],
                    'nr_class_parent_id': parent['name']['class_id']
                })

        return data

    def data(self, pdbs, **kwargs):
        grouping = self.cached('nr')
        if not grouping:
            raise core.InvalidState("No grouping loaded")

        mapping = self.mapping(grouping)
        if not mapping:
            raise core.InvalidState("No mapping loaded")

        return self.parents(grouping, mapping)
