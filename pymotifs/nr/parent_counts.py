"""Import the parent counts.

This will compute the number of changes relative to the last release and put
then in the database.
"""

import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod

from pymotifs.nr.release import Loader as ReleaseLoader
from pymotifs.nr.chains import Loader as ChainLoader
from pymotifs.nr.classes import Loader as ClassLoader
from pymotifs.nr.parents import Loader as ParentLoader


class Loader(core.MassLoader):
    dependencies = set([ReleaseLoader, ChainLoader, ClassLoader, ParentLoader])

    def has_data(self, *args, **kwargs):
        grouping = self.cached('nr')
        if not grouping:
            raise core.Skip("No precomputed grouping to store")

        release_id = grouping[0]['release']
        with self.session() as session:
            query = session.query(mod.NrParentCounts).\
                filter_by(nr_release_id=release_id)
            return bool(query.count())

    def remove(self, *args, **kwargs):
        self.logger.info("No automatic removal of counts or cached data")

    def parent_id(self, nr_id):
        with self.session() as session:
            return session.query(mod.NrReleases).\
                filter_by(nr_release_id=nr_id).\
                one().\
                parent_nr_release_id

    def release_pdbs(self, release_id, resolution):
        with self.session() as session:
            chains = mod.NrChains
            classes = mod.NrClasses
            ife = mod.IfeInfo
            query = session.query(ife.pdb_id).\
                join(chains, chains.ife_id == ife.ife_id).\
                join(classes, classes.nr_class_id == chains.nr_class_id).\
                filter(classes.nr_release_id == release_id).\
                distinct()

            if resolution != 'all':
                query = query.filter(classes.resolution == resolution)

            return [result.pdb_id for result in query]

    def known_names(self, release_id, resolution):
        with self.session() as session:
            classes = mod.NrClasses
            query = session.query(classes.handle, classes.version).\
                filter_by(nr_release_id=release_id, resolution=resolution)
            return [(result.handle, result.version) for result in query]

    def counts(self, parent_id, nr_classes):
        grouper = lambda c: c['name']['cutoff']
        grouped = it.groupby(sorted(nr_classes, key=grouper), grouper)
        release = nr_classes[0]['release']

        data = []
        for resolution, klasses in grouped:
            counts = coll.defaultdict(int)
            parent_names = set(self.known_names(parent_id, resolution))
            current_names = set()
            for klass in klasses:
                counts[klass['name']['type']] += 1
                current_names.add((klass['name']['handle'],
                                   klass['name']['version']))

            parent_pdbs = set(self.release_pdbs(parent_id, resolution))
            current_pdbs = set(self.release_pdbs(release, resolution))

            data.append({
                'resolution': resolution,
                'nr_release_id': release,
                'parent_nr_release_id': parent_id,
                'new_class_count': counts['new'],
                'updated_class_count': counts['updated'],
                'removed_class_count': len(parent_names - current_names),
                'unchanged_class_count': counts['exact'],
                'pdb_added_count': len(current_pdbs - parent_pdbs),
                'pdb_removed_count': len(parent_pdbs - current_pdbs)
            })
        return data

    def data(self, *args, **kwargs):
        grouping = self.cached('nr')
        if not grouping:
            raise core.InvalidState("No grouping loaded")

        parent = self.parent_id(grouping[0]['release'])
        if not parent:
            raise core.Skip("No parent release")

        return [mod.NrParentCounts(**c) for c in self.counts(parent, grouping)]
