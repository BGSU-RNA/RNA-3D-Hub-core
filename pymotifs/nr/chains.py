import datetime as dt

from pymotifs import core

from pymotifs import models as mod

from pymotifs.nr.builder import Builder
from pymotifs.nr.classes import Loader as ClassLoader


class Loader(core.MassLoader):
    dependencies = set([ClassLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days
    table = mod.NrChains

    def has_data(self, *args, **kwargs):
        grouping = self.cached('nr')
        if not grouping:
            raise core.Skip("No precomputed grouping to store")

        release_id = grouping[0]['release']
        with self.session() as session:
            query = session.query(mod.NrChains).\
                join(mod.NrClasses, mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
                filter(mod.NrClasses.nr_release_id == release_id)

            return bool(query.count())

    def remove(self, *args, **kwargs):
        self.logger.info("No automatic removal of chains or cached data")

    def mapping(self, grouping):
        helper = Builder(self.config, self.session)
        release_id = grouping[0]['release']
        classes = [g['name']['full'] for g in grouping]
        return helper.class_id_mapping(classes, release_id)

    def chains(self, grouping, mapping):
        if not grouping:
            raise core.InvalidState("Cannot load chains without classes")

        if not mapping:
            raise core.InvalidState("Cannot load chains without name mapping")

        data = []
        for group in grouping:
            for chain in group['members']:
                if group['name']['full'] not in mapping:
                    raise core.InvalidState("Group %s not in mapping" % group)

                data.append({
                    'ife_id': chain['id'],
                    'nr_class_id': mapping[group['name']['full']],
                    'nr_release_id': group['release'],
                    'rank': chain['rank'],
                    'rep': chain['rank'] == 0
                })

        return data

    def data(self, *args, **kwargs):
        grouping = self.cached('nr')
        if not grouping:
            raise core.InvalidState("No grouping loaded")

        mapping = self.mapping(grouping)
        return self.chains(grouping, mapping)
