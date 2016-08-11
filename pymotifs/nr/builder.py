"""This module contains the code for creating a new NR set. It can find and
group all ifes from structures into classes and then pick a representative of
that class.
"""

import copy
import operator as op
import itertools as it
import collections as coll

from pymotifs import core

from pymotifs.models import NrChains
from pymotifs.models import NrClasses

from pymotifs.utils.naming import Namer
from pymotifs.utils.naming import ChangeCounter
from pymotifs.nr.groups.simplified import Grouper

from pymotifs.constants import NR_BP_PERCENT_INCREASE
from pymotifs.constants import NR_LENGTH_PERCENT_INCREASE
from pymotifs.constants import RESOLUTION_GROUPS
from pymotifs.constants import NR_CLASS_NAME


class Known(core.Base):
    def handles(self):
        """Return a set of all known handles for the nr set.
        """
        with self.session() as session:
            query = session.query(NrClasses.handle).distinct()
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
            query = session.query(NrChains.ife_id.label('id'),
                                  NrClasses.handle.label('handle'),
                                  NrClasses.version.label('version'),
                                  NrClasses.nr_class_id.label('class_id'),
                                  NrClasses.name.label('full_name'),
                                  ).\
                join(NrClasses,
                     NrChains.nr_class_id == NrClasses.nr_class_id).\
                filter(NrClasses.nr_release_id == release_id).\
                filter(NrClasses.resolution == cutoff).\
                order_by(NrClasses.nr_class_id)

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
                mapping[result.name] = result.nr_class_id

        if len(mapping) != len(names):
            raise core.InvalidState("Could not map all names")

        return mapping


class Builder(core.Base):
    """Class to build a new nr set. This handles the logic of grouping the ifes
    into groups, filter them by resolutions, name them using the previous
    release as well as determine the representative.
    """

    def group(self, pdbs, use_discrepancy=True, **kwargs):
        """Group all pdbs into nr sets.
        """
        if not pdbs:
            raise core.InvalidState("Must give pdbs to group")
        grouper = Grouper(self.config, self.session)
        grouper.use_discrepancy = use_discrepancy
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

        self.logger.info("Building nr release with %i pdbs", len(pdbs))

        known = Known(self.config, self.session)
        groups = self.group(pdbs)
        sequence_only = self.group(pdbs, use_discrepancy=False)

        parents = {}
        for cutoff in cutoffs:
            parents[cutoff] = known.classes(parent_release, cutoff)

        named = self.named(groups, parents['all'])
        filtered = self.filter_groups(named, cutoffs)
        with_reps = self.find_representatives(filtered)

        return {
            'parent_counts': self.counts(parents, with_reps),
            'groups': with_reps,
            'sequence_only': sequence_only,
            'release': current_release,
            'parent': parent_release,
        }


class RepresentativeFinder(core.Base):
    """A class to find the representative for a group of ifes. This will find
    the best in terms of bp/nt and attempt to find all those with more bp's and
    nts in the set.
    """

    def sorting_key(self, chain):
        """Function to use for sorting by bps/nt. Deals with things where there
        are 0 bps or nts. It must have bps, length and id entry.

        :param dict chain: Chain to compute a sorting key for.
        :returns: A tuple that can be used to sort the chains.
        """

        ratio = 0
        if chain['bp'] and chain['length']:
            ratio = float(chain['bp']) / float(chain['length'])
        resolution = chain.get('resolution')
        if resolution:
            resolution = resolution * -1
        return (ratio, resolution, chain['id'])

    def naive_best(self, group):
        """Find the best chain terms of bps/nts. This is the starting point for
        finding the representative in a set of ifes. This method is naive
        because it does not favor more complete structures. In addition, it is
        very sensitive to minor changes in number of basepairs and nts.

        :param list group: A list of dictonaries to find the naive
        representative of.
        :returns: The initial representative.
        """
        return max(group, key=self.sorting_key)

    def candidates(self, best, group):
        """Find all possible candidates for a representative within the group,
        given a current best ife. This finds all chains that have at least as
        many basepairs and nucleotides as the best chain. The chains will be
        returned in sorted order.

        :param dict best: The current best.
        :param list group: The list of dicts to search.
        :returns: The list of candidates for the representative.
        """

        len = op.itemgetter('length')
        bp = op.itemgetter('bp')
        same = lambda c: bp(c) == bp(best) and len(c) == len(best)
        possible = it.ifilter(lambda c: len(c) >= len(best), group)
        possible = it.ifilter(lambda c: bp(c) >= bp(best), possible)
        possible = it.ifilterfalse(same, possible)
        return sorted(possible, key=self.sorting_key)

    def increase(self, first, second, key):
        """Compute the percent increase for the given set of dictionaries and
        with the given key. If the second one is 0 then we return 100 for 100%
        increase.

        :param dict first: Dictionary to get the increase to.
        :param dict second: Dictionary to get the increase from.
        :param str key: Key to use
        :returns: The percent increase.
        """

        if not second[key]:
            if not first[key]:
                return 0
            return 100
        return (float(first[key]) / float(second[key]) - 1) * 100

    def best_above_cutoffs(self, representative, candidates):
        """This will find the true representative given a current one and a
        list of candidates. This will attempt to maximize the number of bps and
        nts in the candidate as compared to the current representative. In
        addition, it will only change representatives if we have have enough of
        an increase. This adds stability to the process so minor improvements
        are ignored, while large ones will lead to large changes.

        :param dict representative: The current representative.
        :param list candidates: A list of candidates to examine.
        :returns: The new representative.
        """

        cutoff = (NR_LENGTH_PERCENT_INCREASE, NR_BP_PERCENT_INCREASE)
        possible = []
        for candidate in candidates:
            length_change = self.increase(candidate, representative, 'length')
            bp_change = self.increase(candidate, representative, 'bp')
            if (length_change, bp_change) >= cutoff:
                possible.append(candidate)

        if not possible:
            return representative
        return max(possible, key=self.sorting_key)

    def __call__(self, possible):
        """Find the representative for the group.

        :param list group: List of ifes to find the best for.
        :returns: The ife which should be the representative.
        """

        if not possible:
            raise core.InvalidState("No ifes given")

        # Prefer any xray over any cyro em, as cyro modesl are generally built
        # using x-ray and not yet carefully modeled.
        group = [ife for ife in possible if ife['method'] == 'X-RAY DIFFRACTION']
        if not group:
            group = possible

        best = self.naive_best(group)
        if not best:
            raise core.InvalidState("No current representative")
        self.logger.debug("Naive representative: %s", best['id'])

        candidates = self.candidates(best, group)
        self.logger.debug("Found %i representative candidates",
                          len(candidates))

        rep = self.best_above_cutoffs(best, candidates)
        if not rep:
            raise core.InvalidState("No representative found")

        if rep['id'] != best['id']:
            self.logger.info("Changed representative from %s to %s",
                             best['id'], rep['id'])

        self.logger.debug("Computed representative: %s", rep['id'])
        return rep
