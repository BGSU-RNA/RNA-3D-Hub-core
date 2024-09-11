"""Modify cached data to include class_id.

Other stages in the nr set will need to use the database ids when writing to
the database. Instead of having to keep looking up the mapping and check if
entries exist we put it all in one place.
"""

import copy

from pymotifs import core

from pymotifs.nr.builder import Known
from pymotifs.nr.classes import Loader as ClassLoader
from pymotifs.nr.release import Loader as NrReleaseLoader
from pymotifs.constants import NR_CACHE_NAME


class Loader(core.MassLoader):
    dependencies = set([NrReleaseLoader, ClassLoader])
    allow_no_data = True
    mark = False

    def mapping(self, release_id, grouping):
        """Compute the mapping from names to ids in the database for the given
        nr release id.

        Parameters
        ----------
        release_id : int
            The nr release id
        grouping : list
            The data to comptue a mapping for. Each entry in the list should be
            a dict with a 'name'.'full' entry which contains the NR name.

        Returns
        -------
        mapping : dict
            A mapping from NR name to db id.
        """

        helper = Known(self.config, self.session)
        classes = [g['name']['full'] for g in grouping]
        return helper.mapping(release_id, classes)

    def remove(self, *args, **kwargs):
        """We never remove any data."""
        return False

    def has_data(self, *args, **kwargs):
        """We always run this stage, so this returns False.

        Returns
        -------
        missing : bool
            False
        """
        return False

    def transform(self, grouping, mapping):
        """Transform the grouping to include database ids.

        Parameters
        ----------
        grouping : dict
            The dictonary which must have a 'members' entry to transform
        mapping : dict
            A mapping from NR name to database id.

        Raises
        ------
        InvalidState
            If grouping or mapping is not truthy.

        Returns
        -------
        transformed : dict
            A copy of grouping with database ids in ['name']['class_id'].
        """

        if not grouping:
            raise core.InvalidState("Cannot load chains without classes")

        if not mapping:
            raise core.InvalidState("Cannot load chains without name mapping")

        transformed = copy.deepcopy(grouping)
        for group in transformed:
            for chain in group['members']:
                name = group['name']['full']
                if name not in mapping:
                    raise core.InvalidState("Group %s not in mapping" % group)
                group['name']['class_id'] = mapping[name]

        return transformed

    def data(self, *args, **kwargs):
        """
        Modify the cached data to include database ids. This will load the
        cached NR data and modify it to include the database ids. This will
        never return any data because it is then cached after modification.
        """

        data = self.cached(NR_CACHE_NAME)
        if not data:
            raise core.Skip("Nothing to do here, maybe too few files, maybe an earlier stage failed")

        mapping = self.mapping(data['release'], data['groups'])
        data['groups'] = self.transform(data['groups'], mapping)
        self.cache(NR_CACHE_NAME, data)
        return None
