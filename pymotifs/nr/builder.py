"""
This module contains the code for creating a new NR set. It can find and
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
from pymotifs.nr.groups.simplified import ranking_key

from pymotifs.constants import NR_REPRESENTATIVE_METHOD  # could be "compscore"
from pymotifs.constants import RESOLUTION_GROUPS
from pymotifs.constants import NR_CLASS_NAME
from pymotifs.constants import DNA_CLASS_NAME
from sqlalchemy import desc
from sqlalchemy import func
from pymotifs.constants import WRITE_ALL_EQUIVALENCE_CLASS_RANKINGS


class Known(core.Base):
    def handles(self, molecule_type):
        """Return a set of all known handles for the nr set.
        """
        if molecule_type == 'DNA':
            query_key = molecule_type
        else:
            query_key = 'NR_'
        with self.session() as session:
            query = session.query(mod.NrClasses.handle).\
            filter(func.substr(mod.NrClasses.name, 1, 3) == query_key).\
            distinct()
            return set(result.handle for result in query)

    def classes(self, release_id, cutoff, molecule_type):
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
                'representative': None,
                'name': {
                    'class_id': None,
                    'full': None,
                    'handle': None,
                    'version': None,
                    'cutoff': cutoff,
                },
                'release': release_id,
            }

        def as_member(result):
            return {
                'id': result.id,
                'length': result.length,
                'bp': result.bp,
                'resolution': result.resolution,
            }

        if molecule_type == 'DNA':
            query_key = molecule_type
        else:
            query_key = 'NR_'

        with self.session() as session:
            query = session.query(mod.NrClassRank.ife_id.label('id'),
                                  mod.NrClasses.handle.label('handle'),
                                  mod.NrClasses.version.label('version'),
                                  mod.NrClasses.nr_class_id.label('class_id'),
                                  mod.NrClasses.name.label('full_name'),
                                  mod.NrClassRank.rank,
                                  mod.IfeInfo.bp_count.label('bp'),
                                  mod.IfeInfo.length.label('length'),
                                  mod.PdbInfo.resolution,
                                  ).\
                join(mod.NrClasses,
                     mod.NrClassRank.nr_class_name == mod.NrClasses.name).\
                join(mod.IfeInfo,
                     mod.IfeInfo.ife_id == mod.NrClassRank.ife_id).\
                join(mod.PdbInfo,
                     mod.PdbInfo.pdb_id == mod.IfeInfo.pdb_id).\
                filter(mod.NrClasses.nr_release_id == release_id).\
                filter(mod.NrClasses.resolution == cutoff).\
                filter(mod.IfeInfo.new_style == True).\
                filter(func.substr(mod.NrClasses.name, 1, 3) == query_key).\
                order_by(mod.NrClasses.nr_class_id, mod.NrClassRank.rank)
            self.logger.info("finish the query of the class function in the nr/builder.py")
            results = coll.defaultdict(empty)
            for result in query:
                member = as_member(result)
                results[result.class_id]['members'].append(member)
                # if result.rep:
                #     results[result.class_id]['representative'] = member
                if result.rank == 0:
                    results[result.class_id]['representative'] = member

                results[result.class_id]['name'] = {
                    'class_id': result.class_id,
                    'full': result.full_name,
                    'handle': result.handle,
                    'version': int(result.version),
                    'cutoff': cutoff,
                }

        return results.values()


    def Old_classes(self, release_id, cutoff):
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
                'representative': None,
                'name': {
                    'class_id': None,
                    'full': None,
                    'handle': None,
                    'version': None,
                    'cutoff': cutoff,
                },
                'release': release_id,
            }

        def as_member(result):
            return {
                'id': result.id,
                'length': result.length,
                'bp': result.bp,
                'resolution': result.resolution,
            }

        with self.session() as session:
            query = session.query(mod.NrChains.ife_id.label('id'),
                                mod.NrClasses.handle.label('handle'),
                                mod.NrClasses.version.label('version'),
                                mod.NrClasses.nr_class_id.label('class_id'),
                                mod.NrClasses.name.label('full_name'),
                                mod.NrChains.rep,
                                mod.IfeInfo.bp_count.label('bp'),
                                mod.IfeInfo.length.label('length'),
                                mod.PdbInfo.resolution,
                                ).\
                join(mod.NrClasses,
                    mod.NrChains.nr_class_id == mod.NrClasses.nr_class_id).\
                join(mod.IfeInfo,
                    mod.IfeInfo.ife_id == mod.NrChains.ife_id).\
                join(mod.PdbInfo,
                    mod.PdbInfo.pdb_id == mod.IfeInfo.pdb_id).\
                filter(mod.NrClasses.nr_release_id == release_id).\
                filter(mod.NrClasses.resolution == cutoff).\
                filter(mod.IfeInfo.new_style == True).\
                order_by(mod.NrClasses.nr_class_id, mod.NrChains.rank)

            results = coll.defaultdict(empty)
            for result in query:
                member = as_member(result)
                results[result.class_id]['members'].append(member)
                if result.rep:
                    results[result.class_id]['representative'] = member

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

    def load_parents(self, parent_release, cutoffs, molecule_type):
        parents = {}
        known = Known(self.config, self.session)
        for cutoff in cutoffs:
            parents[cutoff] = known.classes(parent_release, cutoff, molecule_type)
        return parents

    def group(self, pdbs, **kwargs):
        """
        Group all pdbs into equivalence classes.
        This will look at the 'nr' key in the
        config dict to determine if grouping should use discrepancy (key:
        use_discrepancy, default True), use species (key: use_species, default
        True) and enforce single species per group (key: enforce_species,
        default: True).
        This uses the grouper in pymotifs.nr.groups.simplified to group.

        Parameters
        ----------
        pdbs : list
            List of PDB ids to group.

        Returns
        -------
        groups : list
            List of grouped IFE's.
        """

        if len(pdbs) < 100:
            raise core.Skip("Too few pdb files being processed to run the nr stage and make a representative set: %i" % len(pdbs))

        # In July 2024, we only run DNA structures through the command line
        molecule_parent_current = kwargs.get('nr_molecule_parent_current','')
        if molecule_parent_current:
            molecule_type = molecule_parent_current.split(",")[0]
        else:
            molecule_type = 'RNA'


        grouper = Grouper(self.config, self.session)

        if molecule_type == 'RNA':
            self.logger.info('Grouping RNA IFEs into equivalence classes.')
            conf = self.config['nr']
            self.logger.info("conf: %s" % str(conf))
            grouper.use_discrepancy = conf.get('use_discrepancy',
                                            grouper.use_discrepancy)
            grouper.use_species = conf.get('use_species', grouper.use_species)
            grouper.molecule_type = 'RNA'
            if not grouper.use_species:
                enforce = conf.get('enforce_species',
                                grouper.must_enforce_single_species)
                grouper.must_enforce_single_species = enforce
            return grouper(pdbs, **kwargs)
        else:
            self.logger.info('Grouping DNA IFEs into equivalence classes.')
            self.config['dna'] = {'use_discrepancy': True}
            conf = self.config['dna']
            self.logger.info("conf: %s" % str(conf))
            grouper.use_discrepancy = conf.get('use_discrepancy',
                                            grouper.use_discrepancy)
            # grouper.use_species = conf.get('use_species', grouper.use_species)
            # try again to avoid using species
            grouper.molecule_type = 'DNA'
            grouper.use_species = False
            # if not grouper.use_species:
            #     enforce = conf.get('enforce_species',
            #                     grouper.must_enforce_single_species)
            #     ## we do not want to build dna group by single species. Thus, change True to False
            #     ## pipeline have been killed.
            #     ## I think there is too much to group.
            #     ## Thus, Ekko thinks it is better to group by species for test runnings
            #     grouper.must_enforce_single_species = enforce
            return grouper(pdbs, **kwargs)

    def name_groups(self, groups, parents, molecule_type):
        """
        Compute the naming and parents of all groups. Note that this will
        assign parent entries to all groups. However, these parent entries will
        be incorrect when the group is filtered to only those entries within
        a resolution cutoff. This can be corrected by `assign_parents`.

        In other words, each group in groups is a set of ifes from the current release.
        In the previous release, some of these ifes may have been in more than one group.
        Any group an ife was preiously in is called a parent of group.
        Each parent has a handle like 38294 and a version.

        Parameters
        ----------
        groups : list
            List of groups to assign parents to
        parents : list
            List of possible parent groups.

        Returns
        -------
        named : list
            List of groups with 'parents' entries and 'name' filled out.
        """

        ife_count = it.imap(lambda g: len(g['members']), groups)
        ife_count = sum(ife_count)

        self.logger.info("Naming %i groups with %i ifes",
                         len(groups), ife_count)

        namer = Namer(self.config, self.session)
        known = Known(self.config, self.session)
        named = namer(groups, parents, known.handles(molecule_type))

        named_count = it.imap(lambda g: len(g['members']), named)
        named_count = sum(named_count)
        if named_count != ife_count:
            raise core.InvalidState("Missing named ifes")

        self.logger.info("Named %i groups with %i ifes",
                         len(named), named_count)

        return named


    def counts(self, parents, current):
        """
        Compare two releases and count the number of changes.

        This seems to be pretty fragile since any number of things we
        have tried to speed up this code leads to very high change counts.
        """

        data = []
        grouped = it.groupby(current, lambda g: g['name']['cutoff'])
        for cutoff, groups in grouped:
            entry = {'cutoff': cutoff}
            entry.update(self.cutoff_counts(parents[cutoff], list(groups)))
            data.append(entry)
        return data

    def cutoff_counts(self, parents, groups):
        """
        Compute the counts of changes to relative to the parent class. The
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
        """
        Filter the group to produce a new one where all members of the group
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

    def class_name(self, resolution, entry, molecule_type):
        """Create the name of the class given the class and resolution.

        :param dict entry: The class.
        :param str resolution: The resolution cutoff.
        :returns: The class name.
        """
        if molecule_type != 'DNA':
            return NR_CLASS_NAME.format(
                resolution=resolution,
                handle=entry['name']['handle'],
                version=entry['name']['version']
            )
        else:
            return DNA_CLASS_NAME.format(
                resolution=resolution,
                handle=entry['name']['handle'],
                version=entry['name']['version']
            )

    def filter_groups(self, groups, resolutions, molecule_type):
        """This will filter each group so that it contains only the members
        with that pass the given resolution cutoff. The cutoffs should be
        be values that are accepted by within_cutoff. In addition, the groups
        will have their group['name']['full'] and group['name']['cutoff']
        entries set.

        Parameters
        ----------
        groups : list
            The list of group dictonaries to filter.
        resolutions : list
            A list of resolution cutoffs to apply. The cutoffs should be as

        Returns
        -------
        groups : list
            List of groups that have been filtered to contain only the
            requested members.
        """

        # put 'all' first and then decreasing order of resolution; might not work
        # resolutions = sorted(resolutions, key=lambda x: float(x) if x != 'all' else 999.0, reverse=True)

        data = []
        # sort group by handle so they clearly run from 00000 to 99999; might not work
        # for entry in sorted(copy.deepcopy(groups), key=lambda x: x['name']['handle']):
        for entry in copy.deepcopy(groups):
            for resolution in resolutions:
                filtered = self.within_cutoff(entry, resolution)
                if not filtered:
                    continue
                filtered['name']['full'] = self.class_name(resolution, entry, molecule_type)
                filtered['name']['cutoff'] = resolution
                data.append(filtered)

        # it seems to be very important to keep groups together by cutoff
        # maybe because later they will be grouped by cutoff when being recorded
        # still, we can sort from "all" down to 1.5
        # that way, we can use the ranking in the "all" category to determine the
        # ranking in the categories with fewer structures
        # try to also sort by handle ... does not seem to work! Causes all new EC names.
        return sorted(data, key=lambda g: g['name']['cutoff'], reverse=True)


    def attach_parents(self, groups, parents):
        """
        Load all parents of the given representatives. These groups should
        already be named and thus have a parent assignment. However, we need to
        load the parent information for all resolution cutoffs, as well as load
        the representative of the parents, which are not already loaded.

        Parameters
        ----------
        groups : list
            List of all groups to attach the parents to.
        parents : list
            List of all known parents

        Returns
        -------
        groups : list
            List of groups with the 'parents' entries set correctly.
        """

        # return a tuple of cutoff like 3.5 and group handle like 31482
        def as_key(group, cutoff=None):
            if cutoff:
                return (cutoff, group['name']['handle'])
            return (group['name']['cutoff'], group['name']['handle'])

        flattened = it.chain.from_iterable(parents.values())

        # map (cutoff, handle) to parent group, which is a dictionary
        mapping = {as_key(p): p for p in flattened}

        data = []
        for group in copy.deepcopy(groups):
            cutoff = group['name']['cutoff']
            known_parents = []
            for parent in group['parents']:
                key = as_key(parent, cutoff=cutoff)
                if key in mapping:
                    known_parents.append(mapping[key])
            group['parents'] = known_parents
            data.append(group)
        return data


    def using_old_rank(self, group):
        ife_id_list = []
        for parent in group['parents']:
            for member in parent['members']:
                ife_id_list.append(member['id'])
        with self.session() as session:
            query = session.query(mod.NrCqs.ife_id,mod.NrCqs.nr_name,mod.NrCqs.composite_quality_score.label('cqs')).\
            filter(mod.NrCqs.ife_id.in_(ife_id_list)).\
                order_by(desc(mod.NrCqs.nr_name)).distinct(mod.NrCqs.ife_id)
        last_cqs_values = {}
        for result in query:
            last_cqs_values[result.ife_id] = result.cqs

        # self.logger.info("Found last_cqs_values %s", last_cqs_values)

        sorted_last_cqs_values = sorted(last_cqs_values.items(), key=lambda x: x[1])

        # self.logger.info("Found sorted_last_cqs_values %s", sorted_last_cqs_values)

        old_rank = {ife_id: rank for rank, (ife_id, _) in enumerate(sorted_last_cqs_values)}
        ### will return a dict, key is the ife_id and value is the rank of cqs {'ife_id1':0,'ife_id2':1,.....}
        return old_rank

    def find_representatives(self, groups, molecule_type, sorting_key=ranking_key):
        """Compute the representative for each group. This will modify the
        group to now have a 'representative' entry containing the
        representative entry. In addition, the members will be sorted and
        assigned a rank based upon the sorting. The representative will always
        be placed in the first position and thus have rank 0, the lowest
        rank.

        Parameters
        ----------
        groups : list
            List of groups to get representatives for
        sorting_key : function, default ranking_key
            A function to use when sorting the members

        Returns
        -------
        groups : list
            The list of groups which have been modified to include the
            representative entry and the members are resorted.
        """
        data = []
        rep_finder = RepresentativeFinder(self.config, self.session)

        if molecule_type == 'DNA':
            query_key = molecule_type
        else:
            query_key = 'NR_'

        with self.session() as session:
            query = session.query(mod.NrClassRank.nr_class_name).\
            filter(func.substr(mod.NrClassRank.nr_class_name, 1, 3) == query_key).\
            distinct()
        nr_class_name_list = [row.nr_class_name for row in query]

        handle_to_ordered_members = {}
        # self.logger.info("Found groups[0] %s", groups[0])
        # must process groups in the same order as given, apparently order is important
        for group in copy.deepcopy(groups):
            handle = group['name']['handle']

            # in order from all, 4.0, 3.5, 3.0, 20.0, etc.

            # old_rank = self.using_old_rank(group)
            # self.logger.info("Found old_rank %s", old_rank)
            # self.logger.info("Found group %s", group)
            # self.logger.info("Found group_parents_members %s", group['parents'][0]['members'])
            # ## checking if we have new ife for this group
            # if len(old_rank) == len(group['parents'][0]['members']):
            #     # Sort the members based on the rank value in old_rank
            #     sorted_parent_members = sorted(group['parents'][0]['members'], key=lambda x: old_rank[x['id']])
            #     sorted_members = sorted(group['members'], key=lambda x: old_rank[x['id']])
            #     group['members'] = sorted_members
            #     group['representative'] = sorted_members[0]
            #     group['parents'][0]['members'] = sorted_parent_members
            # else:
            ## new if statement
            ## make two log info messages and show what is old and what is new.
            # if not WRITE_ALL_EQUIVALENCE_CLASS_RANKINGS and group['comment'] == 'Exact match':
            if (group['name']['full'] in nr_class_name_list) and (not WRITE_ALL_EQUIVALENCE_CLASS_RANKINGS):
                # self.logger.info("Already have ranking and representative for %s", group['name']['full'])
                # The lines below signal that we should retrieve the previous ranking and rep
                # from the database
                group['representative'] = None
                group['members'] = []
            elif handle in handle_to_ordered_members:
                # Use the ranking from equivalence classes at a higher resolution cutoff
                # Saves time because you don't have to compute CQS all over again
                self.logger.info("Using new ranking for %s", group['name']['full'])
                keep_members = []
                for member in handle_to_ordered_members[handle]:
                    # if mem has resolution less than equal to the current cutoff, keep it
                    # convert the resolution to float in a way that will work even if value is 'all'
                    try:
                        res = float(member['resolution'])
                    except:
                        res = 999999.0
                    if res <= float(group['name']['cutoff']):
                        # use deepcopy to avoid changing the member from previous cutoff
                        keep_members.append(copy.deepcopy(member))

                group['representative'] = keep_members[0]

                self.logger.info('Representative chain %s' % str(group['representative']['id']))

                group['members'] = keep_members
                for index, member in enumerate(group['members']):
                    member['rank'] = index
                    self.logger.info("rank %3d ife %s" % (index, member['id']))
            else:
                # Look up data, calculate CQS, and rank the members
                #
                self.logger.info("Ranking and selecting representative for %s", group['name']['full'])
                ordered_members = rep_finder(group)
                group['representative'] = ordered_members[0]

                self.logger.info('Representative chain %s' % str(group['representative']['id']))

                group['members'] = ordered_members
                for index, member in enumerate(group['members']):
                    member['rank'] = index

                # uncomment the next line to use the elif statement above and save some time
                # but maybe this scrambles the connection to the parents and thus backfires
                handle_to_ordered_members[handle] = ordered_members
            data.append(group)
        return data


    def __call__(self, pdbs, parent_release, current_release,
                 cutoffs=RESOLUTION_GROUPS, **kwargs):
        """
        Build equivalence classes and the representative set.

        :pdbs: The list of pdbs to process.
        :current: Current release id.
        :new: Id for the next release.
        :resolutions: Resolution groups to create a class for.
        :returns: A list of nr classes with their memebers and parents.
        """

        if not pdbs:
            raise core.InvalidState("Must give pdbs to group")

        self.logger.info("Building nr release with %d pdbs" % len(pdbs))

        groups = self.group(pdbs, **kwargs)

        # In July 2024, we only run DNA structures through the command line
        molecule_parent_current = kwargs.get('nr_molecule_parent_current','')
        if molecule_parent_current:
            molecule_type = molecule_parent_current.split(",")[0]
        else:
            molecule_type = 'RNA'

        parents = self.load_parents(parent_release, cutoffs, molecule_type)
        named = self.name_groups(groups, parents['all'], molecule_type)
        filtered = self.filter_groups(named, cutoffs, molecule_type)
        with_parents = self.attach_parents(filtered, parents)

        # self.logger.info("Found with_parents %s", with_parents)
        with_reps = self.find_representatives(with_parents, molecule_type)
        # self.logger.info("Found find_representatives %s", with_reps)

        return {
            'parent_counts': self.counts(parents, with_reps),
            'groups': with_reps,
            'release': current_release,
            'parent': parent_release,
        }


class RepresentativeFinder(core.Base):
    """
    A class to find the representative for a group of ifes.
    This will find the best in terms of the criterion
    NR_REPRESENTATIVE_METHOD.

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
        """Get the method with the given name.

        Parameters
        ----------
        name : str
            The name of the method

        Returns
        -------
            An object that can be called to find the representative.
        """
        # self.logger.info("tracking numbers #0001")
        finder = reps.fetch(name)
        # self.logger.info("tracking numbers #0002")
        ## Actually, the return info is equal to CompSocre(self.config, self.session)
        ## We normally rename a class or initial a class in this way, finder = CompScore(self.config, self.session)
        return finder(self.config, self.session)


    def __call__(self, group, method=NR_REPRESENTATIVE_METHOD):
        """Find the representative for the group.

        Parameters
        ----------
        group : dict
            The grouping to find a representative for.
        method : str, default `pymotifs.constants.NR_REPRESENTATIVE_METHOD`
            Name of the method to use for finding a representative.

        Returns
        -------
            The ife which should be the representative.
        """

        if not group:
            raise core.InvalidState("No group given")

        if method not in self.methods:
            raise core.InvalidState("Unknown method %s" % method)

        finder = self.method(method)
        # self.logger.info("Found group in RepresentativeFinder class: %s", group)
        new_data=finder(group)
        # self.logger.info("Found finder(group) in RepresentativeFinder class: %s", new_data)
        return new_data
