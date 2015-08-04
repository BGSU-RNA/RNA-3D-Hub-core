import datetime as dt

from pymotifs import core

from pymotifs.models import NrChains
from pymotifs.models import NrClasses

from pymotifs.utils import tmp
from pymotifs.nr.classes import Loader as ClassLoader


class Loader(core.MassLoader):
    dependencies = set([ClassLoader])
    update_gap = dt.timedelta(7)  # Only update every 7 days

    def ids(self, names, release_id):
        with self.session() as session:
            query = session.query(NrClasses.id,
                                  NrClasses.name,
                                  ).\
                filter(NrClasses.name.in_(names)).\
                filter(NrClasses.nr_release_id == release_id)

            if query.count() == 0:
                raise core.InvalidState("Found no clases with given names")

            mapping = {}
            for result in query:
                mapping[result.name] = result.id

        if len(mapping) != len(names):
            raise core.InvalidState("Could not find all names")

        return mapping

    def chains(self, grouping, mapping):
        if not grouping:
            raise core.InvalidState("Cannot load chains without classes")

        if not mapping:
            raise core.InvalidState("Cannot load chains without mapping")

        data = []
        for (name, release_id), chains in grouping.items():

            if name not in mapping:
                raise core.InvalidState("Could not map name %s" % name)

            class_id = mapping[name]
            for chain in chains:
                data.append({
                    'autonomous_group_id': chain['id'],
                    'nr_class_id': class_id,
                    'nr_release_id': release_id,
                    'rank': chain['rank'],
                    'rep': chain['rank'] == 0
                })

        return data

    def data(self, *args, **kwargs):
        grouping = tmp.load('nr')
        names = [key[0] for key in grouping.keys()]
        release_id = grouping.keys()[0][1]
        mapping = self.ids(names, release_id)

        return [NrChains(**chain) for chain in self.chains(grouping, mapping)]
