"""This module contains the code for creating a new NR set. It can find and
group all ifes from structures into classes and then pick a representative of
that class.
"""

import copy
import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod

from pymotifs.utils.naming import Namer
from pymotifs.nr import representatives as reps
from pymotifs.utils.naming import ChangeCounter
from pymotifs.nr.groups.simplified import Grouper

from pymotifs.constants import NR_REPRESENTATIVE_METHOD
from pymotifs.constants import RESOLUTION_GROUPS
from pymotifs.constants import NR_CLASS_NAME


class Known(core.Base):
    def handles(self):
        """Return a set of all known handles for the nr set.
        """
        with self.session() as session:
            query = session.query(mod.NrClasses.handle).distinct()
            return set(result.handle for result in query)

    def classes(self, release_id, cutoff):
        """Get all classes with the given resolution cutoff in the given
        release. If nothing is below the cutoff then an empty dictonary will be
        returned.

        :release_id: Release id to lookup.
        :cutoff: The resolution cutoff to use. Default is 'all'.
        :returns A list of dictonaries for each class.
        """

        def empty():
            return {
                'members': [],
                'name': {
                    'class_id': None,
                    'full': None,
                    'handle': None,
                    'version': None,
                    'cutoff': cutoff,
                },
                'release': release_id,
            }

        with self.session() as session:
            query = session.query(mod.NrChains.ife_id.label('id'),
                                  mod.NrClasses.handle.label('handle'),
                                  mod.NrClasses.version.label('version'),
                                  mod.NrClasses.nr_class_id.label('class_id'),
                                  mod.NrClasses.name.label('full_name'),
                                  ).\
                join(mod.NrClasses,
                     mod.NrChains.nr_class_id == mod.NrClasses.nr_class_id).\
                filter(mod.NrClasses.nr_release_id == release_id).\
                filter(mod.NrClasses.resolution == cutoff).\
                order_by(mod.NrClasses.nr_class_id)

            results = coll.defaultdict(empty)
            for result in query:
                results[result.class_id]['members'].append({'id': result.id})
                results[result.class_id]['name'] = {
                    'class_id': result.class_id,
                    'full': result.full_name,
                    'handle': result.handle,
                    'version': int(result.version),
                    'cutoff': cutoff,
                }

        return results.values()

    def mapping(self, release_id, names):
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
            query = session.query(mod.NrClasses.nr_class_id,
                                  mod.NrClasses.name,
                                  ).\
                filter(mod.NrClasses.name.in_(names)).\
                filter(mod.NrClasses.nr_release_id == release_id)

            if query.count() == 0:
                self.logger.info(names)
                raise core.InvalidState("Found no clases with given names")

            mapping = {}
            for result in query:
                mapping[result.name] = result.nr_class_id

        if len(mapping) != len(names):
            raise core.InvalidState("Could not map all names")

        return mapping


class Builder(core.Base):
    """Class to build a new nr set. This handles the logic of grouping the ifes
    into groups, filter them by resolutions, name them using the previous
    release as well as determine the representative.
    """

    def group(self, pdbs, **kwargs):
        """Group all pdbs into nr sets. This will look at the 'nr key in the
        config dict to determine if grouping should use discrepancy (key:
        use_discrepancy, deafult True), use species (key: use_species, default
        True) and enforce single species per group (key: enforce_species,
        default: True). This uses the grouper in pymotifs.nr.groups.simplified
        to group.

        Parameters
        ----------
        pdbs : list
            List of PDB ids to group.

        Returns
        -------
        groups : list
            List of grouped IFE's.
        """
        grouper = Grouper(self.config, self.session)
        conf = self.config['nr']
        grouper.use_discrepancy = conf.get('use_discrepancy',
                                           grouper.use_discrepancy)
        grouper.use_species = conf.get('use_species', grouper.use_species)
        if not grouper.use_species:
            enforce = conf.get('enforce_species',
                               grouper.must_enforce_single_species)
            grouper.must_enforce_single_species = enforce
        return grouper(pdbs, **kwargs)

    def named(self, groups, parents):
        """Compute the naming and parents of all groups.
        """

        ife_count = it.imap(lambda g: len(g['members']), groups)
        ife_count = sum(ife_count)

        self.logger.info("Naming %i groups with %i ifes",
                         len(groups), ife_count)

        namer = Namer(self.config, self.session)
        known = Known(self.config, self.session)
        named = namer(groups, parents, known.handles())

        named_count = it.imap(lambda g: len(g['members']), named)
        named_count = sum(named_count)
        if named_count != ife_count:
            raise core.InvalidState("Missing named ifes")

        self.logger.info("Named %i groups with %i ifes",
                         len(named), named_count)

        return named

    def counts(self, parents, current):
        """Compare two releases and count the number of changes.
        """

        data = []
        grouped = it.groupby(current, lambda g: g['name']['cutoff'])
        for cutoff, groups in grouped:
            entry = {'cutoff': cutoff}
            entry.update(self.cutoff_counts(parents[cutoff], list(groups)))
            data.append(entry)
        return data

    def cutoff_counts(self, parents, groups):
        """Compute the counts of changes to relative to the parent class. The
        groups should be from the named method, while parents should be as from
        Known.motifs.

        :param list groups: The list of groups.
        :param list parents: The parent groups.
        :returns: A dictonary with 'pdbs' and 'classes' entries, which
        summarize the changes for the release. The classes entry the same as
        the 'groups' entry produced by ChangeCounter. And the pdbs entry is the
        result of a transform.
        """

        def as_pdbs(group):
            return [m['id'].split('|')[0] for m in group['members']]

        counter = ChangeCounter(self.config, self.session)
        counts = counter(groups, parents, pdbs=as_pdbs)
        counts['ifes'] = counts.pop('members')
        counts['pdbs'].pop('unchanged')
        counts['classes'] = counts.pop('groups')
        return counts

    def within_cutoff(self, group, cutoff):
        """Filter the group to produce a new one where all members of the group
        are within the resolution cutoff. This will update the rank of the
        members so that they indicate the new rank.

        :param dict group: The group to filter.
        :param str cutoff: The resolution cutoff to use.
        """

        if cutoff == 'all':
            return copy.deepcopy(group)

        cutoff = float(cutoff)
        updated = copy.deepcopy(group)

        filtered = updated['members']
        filtered = it.ifilter(lambda c: c['resolution'] is not None, filtered)
        filtered = it.ifilter(lambda c: c['resolution'] <= cutoff, filtered)
        filtered = list(filtered)
        if not filtered:
            return {}

        for index, entry in enumerate(filtered):
            entry['rank'] = index

        updated['members'] = filtered
        return updated

    def class_name(self, resolution, entry):
        """Create the name of the class given the class and resolution.

        :param dict entry: The class.
        :param str resolution: The resolution cutoff.
        :returns: The class name.
        """

        return NR_CLASS_NAME.format(
            resolution=resolution,
            handle=entry['name']['handle'],
            version=entry['name']['version']
        )

    def filter_groups(self, groups, resolutions):
        data = []
        for entry in groups:
            for resolution in resolutions:
                filtered = self.within_cutoff(entry, resolution)
                if not filtered:
                    continue
                filtered['name']['full'] = self.class_name(resolution, entry)
                filtered['name']['cutoff'] = resolution
                data.append(filtered)
        return sorted(data, key=lambda g: g['name']['cutoff'])

    def find_representatives(self, groups):
        data = []
        rep_finder = RepresentativeFinder(self.config, self.session)
        for group in groups:
            rep = rep_finder(group['members'])
            rep['rank'] = 0
            group['representative'] = rep
            group['members'].remove(rep)
            group['members'].insert(0, rep)
            for index, member in enumerate(group['members']):
                member['rank'] = index
            data.append(group)
        return data

    def __call__(self, pdbs, parent_release, current_release,
                 cutoffs=RESOLUTION_GROUPS, **kwargs):
        """Build the nr set.

        :pdbs: The list of pdbs to process.
        :current: Current release id.
        :new: Id for the next release.
        :resolutions: Resolution groups to create a class for.
        :returns: A list of nr classes with their memebers and parents.
        """

        if not pdbs:
            raise core.InvalidState("Must give pdbs to group")

        self.logger.info("Building nr release with %i pdbs", len(pdbs))

        groups = self.group(pdbs)

        parents = {}
        known = Known(self.config, self.session)
        for cutoff in cutoffs:
            parents[cutoff] = known.classes(parent_release, cutoff)

        named = self.named(groups, parents['all'])
        filtered = self.filter_groups(named, cutoffs)
        with_reps = self.find_representatives(filtered)

        return {
            'parent_counts': self.counts(parents, with_reps),
            'groups': with_reps,
            'release': current_release,
            'parent': parent_release,
        }


class RepresentativeFinder(core.Base):
    """A class to find the representative for a group of ifes. This will find
    the best in terms of bp/nt and attempt to find all those with more bp's and
    nts in the set.

    Attributes
    ----------
    methods : set
        A set of the known method names.
    """

    @property
    def methods(self):
        """Get the known methods for selecting a representative.
        """
        if not hasattr(self, '_methods'):
            self._methods = set(n for n, k in reps.known())
        return self._methods

    def method(self, name):
        """Get a method
        """
        finder = reps.fetch(name)
        return finder(self.config, self.session)

    def __call__(self, possible, method=NR_REPRESENTATIVE_METHOD):
        """Find the representative for the group.

        :param list group: List of ifes to find the best for.
        :returns: The ife which should be the representative.
        """

        if not possible:
            raise core.InvalidState("No ifes given")

        if method not in self.methods:
            raise core.InvalidState("Unknown method %s" % method)

        finder = self.method(method)
        return finder(possible)
