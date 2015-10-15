import datetime as dt
from pprint import pprint

from pymotifs import core

from pymotifs.models import NrChains
from pymotifs.models import NrClasses

from pymotifs.utils import tmp
from pymotifs.nr.builder import Builder
from pymotifs.nr.classes import Loader as ClassLoader


class Loader(core.MassLoader):
    dependencies = set([ClassLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days

    def has_data(self, *args, **kwargs):
        grouping = tmp.load('nr')
        if not grouping:
            raise core.Skip("No precomputed grouping to store")

        release_id = grouping[0]['release']
        with self.session() as session:
            query = session.query(NrChains).\
                join(NrClasses, NrClasses.id == NrChains.nr_class_id).\
                filter(NrClasses.nr_release_id == release_id)

            return bool(query.count())

    def remove(self, *args, **kwargs):
        tmp.cleanup('nr')

    def mapping(self, grouping):
        helper = Builder(self.config, self.session)
        release_id = grouping[0]['release']
        classes = [g['name']['full'] for g in grouping]
        return helper.class_id_mapping(classes, release_id)

    def chains(self, grouping):
        if not grouping:
            raise core.InvalidState("Cannot load chains without classes")

        mapping = self.mapping(grouping)
        if not grouping:
            raise core.InvalidState("Cannot load chains without name mapping")

        data = []
        for group in grouping:
            for chain in group['members']:
                data.append({
                    'ife_group_id': chain['id'],
                    'nr_class_id': mapping[group['name']['full']],
                    'nr_release_id': group['release'],
                    'rank': chain['rank'],
                    'rep': chain['rank'] == 0
                })

        return data

    def data(self, *args, **kwargs):
        grouping = tmp.load('nr')
        return [NrChains(**chain) for chain in self.chains(grouping)]
