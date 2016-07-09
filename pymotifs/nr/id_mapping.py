"""Modify cached data to include class_id.

Other stages in the nr set will need to use the database ids when writing to
the database. Instead of having to keep looking up the mapping and check if
entries exist we put it all in one place.
"""

from pymotifs import core

from pymotifs.nr.builder import Known
from pymotifs.nr.classes import Loader as ClassLoader
from pymotifs.nr.release import Loader as NrReleaseLoader
from pymotifs.constants import NR_CACHE_NAME


class Loader(core.MassLoader):
    dependencies = set([NrReleaseLoader, ClassLoader])
    allow_no_data = True

    def mapping(self, release_id, grouping):
        helper = Known(self.config, self.session)
        classes = [g['name']['full'] for g in grouping]
        return helper.mapping(release_id, classes)

    def remove(self, *args, **kwargs):
        return False

    def has_data(self, *args, **kwargs):
        return False

    def transform(self, grouping, mapping):
        if not grouping:
            raise core.InvalidState("Cannot load chains without classes")

        if not mapping:
            raise core.InvalidState("Cannot load chains without name mapping")

        for group in grouping:
            for chain in group['members']:
                name = group['name']['full']
                if name not in mapping:
                    raise core.InvalidState("Group %s not in mapping" % group)
                group['name']['class_id'] = mapping[name]

        return grouping

    def data(self, *args, **kwargs):
        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.InvalidState("No grouping loaded")

        mapping = self.mapping(data['release'], data['groups'])
        data['groups'] = self.transform(data['groups'], mapping)
        self.cache(NR_CACHE_NAME, data)
        return None
