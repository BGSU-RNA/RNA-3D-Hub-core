import copy
import itertools as it

from pymotifs import core

from pymotifs.models import NrChains
from pymotifs.models import NrClasses

from pymotifs.nr.groups.naming import Namer
from pymotifs.nr.groups.simplified import Grouper

RESOLUTION_GROUPS = ['1.5', '2.0', '2.5', '3.0', '3.5', '4.0', '20.0',
                     'all']


class Builder(core.Base):

    def class_id_mapping(self, names, release_id):
        """Create a mapping from nr class names to id in the database. This
        will raise an exception if it cannot find all names or if not given a
        list of names and a release id.

        :names: A list of names.
        :release_id: The release id to use.
        :returns: A dictonary mapping class_name to id.
        """

        if not names or not release_id:
            raise core.InvalidState("Must give names and release id")

        with self.session() as session:
            query = session.query(NrClasses.nr_class_id,
                                  NrClasses.name,
                                  ).\
                filter(NrClasses.name.in_(names)).\
                filter(NrClasses.nr_release_id == release_id)

            if query.count() == 0:
                self.logger.info(names)
                raise core.InvalidState("Found no clases with given names")

            mapping = {}
            for result in query:
                mapping[result.name] = result.id

        if len(mapping) != len(names):
            raise core.InvalidState("Could not map all names")

        return mapping

    def known_handles(self):
        """Return a set of all known handles for the nr set.
        """

        with self.session() as session:
            query = session.query(NrClasses.handle).distinct()
            return set(result.handle for result in query)

    def group(self, pdbs, **kwargs):
        """Group all pdbs into nr sets.
        """

        grouper = Grouper(self.config, self.session)
        return grouper(pdbs, **kwargs)

    def named(self, groups, parents):
        """Compute the naming and parents of all groups.
        """

        namer = Namer(self.config, self.session)
        handles = self.known_handles()
        return namer(groups, parents, handles)

    def within_cutoff(self, group, cutoff):
        """Filter the group to produce a new one where all members of the group
        are within the resolution cutoff. This will update the rank of the
        members so that they indicate the new rank.

        :param dict group: The group to filter.
        :param str cutoff: The resolution cutoff to use.
        """

        updated = copy.deepcopy(group)
        chains = updated['members']

        filtered = chains
        if cutoff != 'all':
            filtered = it.ifilter(lambda c: c['resolution'] is not None)
            filtered = it.ifilter(lambda c: c['resolution'] <= float(cutoff))
        filtered = list(filtered)

        if not filtered:
            return {}

        filtered.sort(key=lambda f: f['rank'])
        for index, entry in enumerate(filtered):
            entry['rank'] = index

        updated['members'] = filtered
        return updated

    def load_classes_in_release(self, release_id, cutoff='all'):
        """Get all classes with the given resolution cutoff in the given
        release. If nothing is below the cutoff then an empty dictonary will be
        returned.

        :release_id: Release id to lookup.
        :cutoff: The resolution cutoff to use. Default is 'all'.
        :returns A list of dictonaries for each class.
        """

        with self.session() as session:
            query = session.query(NrChains.ife_id.label('id'),
                                  NrClasses.handle.label('handle'),
                                  NrClasses.version.label('version'),
                                  NrClasses.nr_class_id.label('class_id'),
                                  NrClasses.name.label('full_name'),
                                  ).\
                join(NrClasses, NrChains.nr_class_id == NrClasses.nr_class_id).\
                filter(NrClasses.nr_release_id == release_id).\
                filter(NrClasses.resolution == cutoff).\
                order_by(NrClasses.nr_class_id)

            grouped = it.groupby(query, lambda v: (v.handle, v.version))

            results = []
            for (handle, version), members in grouped:
                members = list(members)
                class_id = members[0].class_id
                full_name = members[0].full_name
                results.append({
                    'members': [{'id': member.id} for member in members],
                    'name': {
                        'class_id': class_id,
                        'full': full_name,
                        'handle': handle,
                        'version': int(version)
                    }
                })

        return results

    def class_name(self, entry, resolution):
        return 'NR_{resolution}_{handle}.{version}'.format(
            resolution=resolution,
            handle=entry['name']['handle'],
            version=entry['name']['version']
        )

    def __call__(self, pdbs, current, new, resolutions=RESOLUTION_GROUPS,
                 **kwargs):

        """Build the nr set

        :pdbs: The list of pdbs to process.
        :current: Current release id.
        :new: Id for the next release.
        :resolutions: Resolution groups to create a class for.
        :returns: A list of nr classes with their memebers and parents.
        """

        known = self.load_classes_in_release(current)
        groups = self.group(pdbs)
        named = self.named(groups, known)

        data = []
        for entry in named:
            for resolution in resolutions:
                filtered = self.within_cutoff(entry, resolution)
                if not filtered:
                    continue

                filtered['release'] = new
                filtered['name']['full'] = self.class_name(entry, resolution)
                filtered['name']['cutoff'] = resolution
                data.append(filtered)

        return data
