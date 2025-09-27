"""
Create a new motif release and cluster all motifs. This will
cluster loops from all structure IFE's, that are representatives of their NR
class and come from X-ray or EM structures. The clustered motifs are not yet stored
in the database but are just computed and temporary files are written.

Manual Options
--------------
nr_release : str
    The NR release to use, instead of the newest one.
loop_size_limit : int
    The maximum loop size to use
loop_exclude_file : str
    Name of a file to that contains loops to exclude from clustering. This is
    meant for testing out different excluding procedures without update loop
    quality. The file should contain a single loop id per line, and only that.
loop_type : str
    The specific loop type to process
"""

import collections
import os
from sqlalchemy import desc
from sqlalchemy import or_, and_

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import releases as rel

from pymotifs.constants import MOTIF_ALLOWED_METHODS
from pymotifs.constants import MOTIF_RESOLUTION_CUTOFF

from pymotifs.motif_atlas.builder import Builder

from pymotifs.nr.release import Loader as NrReleaseLoader
# from pymotifs.loops.release import Loader as LoopReleaseLoader

from pymotifs.nr.loader import Loader as NrLoader
from pymotifs.loops.loader import Loader as LoopLoader

Releases = collections.namedtuple('Releases', [
    'nr',
    'parent',
    'parent_index',
    'current'
])

""" A list of loop ids we should not use for clustering. """
BLACKLIST = set([
    'HL_3ICQ_004',
    'HL_3V2F_005',
    'HL_1Y0Q_002',
    'HL_2IL9_002',
    'HL_2IL9_005',
    'IL_4WF9_111'     # 2025-03-27 much too large, near cWW pairs would split it into a J3
])


class Loader(core.MassLoader):
    """
    A loader to create a motif release as well as process and cache all release data.
    """

    types = ['J3', 'IL', 'HL']
    types = ['HL', 'IL', 'J3']
    types = ['IL', 'J3', 'HL']
    types = ['HL', 'IL', 'J3','J4','J5','J6','J7','J8','J9']

    dependencies = set([NrLoader, LoopLoader])
    mark = False       # only works when to_process returns a list of pdb ids
    use_marks = False  # no need to look for marks of what was already done

    merge_data = True  # allow overwriting previously written release data

    def nr_release_id(self, before_date=None, **kwargs):
        """
        Get the Representative Set release number, if not given manually.
        If no before_date is given then we get the latest,
        otherwise we get the release for the given date.
        If no release exists for that date then we fail.

        :param date before_date: The date to use.
        :returns: The nr release id.
        """

        if 'nr_release' in kwargs.get('manual', {}):
            return kwargs['manual']['nr_release']

        if before_date is None:
            nr_release, _ = NrReleaseLoader(self.config, self.session).\
                current_id()
            return nr_release

        with self.session() as session:
            query = session.query(mod.NrReleases).\
                filter_by(date=before_date)

            if query.count() != 1:
                raise core.InvalidState("No nr release on %s", before_date)
            return query.one().nr_release_id


    def get_nr_release_date(self, nr_release_id):
        """
        Get the date of the specified nr release
        """

        with self.session() as session:
            query = session.query(mod.NrReleases).\
                filter_by(nr_release_id=nr_release_id)

            if query.count() != 1:
                raise core.InvalidState("Unclear date for nr release %s", nr_release_id)

            return query.one().date


    def current_id(self, **kwargs):
        """
        Get the current motif release id and the index.
        That is supposed to be a value of ml_release_id ml_release_id that is
        already in the database.
        If there is no release then the release_id is 0.0 and the index is 0.

        :returns: The current release id and index in the ml_releases table
        """

        if 'ml_release_id' in kwargs.get('manual', {}):
            manual_current_id, next_id = kwargs['manual']['ml_release_id'].split(",")
            self.logger.info("We manually get the current ml_release_id: %s" % manual_current_id)
            if manual_current_id == '0.0':
                return ('0.0', 0)
            else:
                try:
                    with self.session() as session:
                        query = session.query(mod.MlReleases.ml_release_id,
                                    mod.MlReleases.index).\
                        order_by(desc(mod.MlReleases.index)).\
                        filter(mod.MlReleases.ml_release_id == manual_current_id).\
                        filter(mod.MlReleases.type == 'HL').\
                        limit(1)

                        current = query.one()
                        self.logger.info('Manually get the ml release id %s' % (manual_current_id))
                        self.logger.info('current from manually defined ml_release_id %s and its index %s' % (current.ml_release_id , current.index))
                        return (current.ml_release_id, current.index)
                except Exception:
                    return (manual_current_id, 0)
        else:
            # get the most recent ml_release_id
            with self.session() as session:
                query = session.query(mod.MlReleases.ml_release_id,
                                    mod.MlReleases.index).\
                    order_by(desc(mod.MlReleases.index)).\
                    filter(mod.MlReleases.type == self.types[0]).\
                    limit(1)

                if query.count() == 0:
                    return ('0.0', 0)

                current = query.one()
                self.logger.info('Current ml_release_id %s and index %s' % (current.ml_release_id , current.index))
                return (current.ml_release_id, current.index)


    def next_motif_release_id(self, **kwargs):
        """
        Get the next value for ml_release_id, the one we are creating now.

        :returns: The next motif id.
        """
        if 'ml_release_id' in kwargs.get('manual', {}):
            # set the next release id manually
            ml_release_ids = kwargs['manual']['ml_release_id'].split(',')
            return ml_release_ids[1]
        else:
            # get the most recent ml_release_id
            motif_release, index = self.current_id()
            # look up whether we are incrementing the major or minor index; like 3.90 to 3.91
            mode = self.config['release_mode']['motifs']
            return rel.next_id(motif_release, mode=mode)


    def to_process(self, pdbs, **kwargs):
        """
        Determine what to process.
        This stage works on a set of nr releases.
        Ignores pdbs, because we pass in nr_release_id manually or look it up.

        Returns
        -------
        releases : Release
            A Release object describing all releases to use.
        """

        # get the most recent existing motif atlas release id
        parent, parent_index = self.current_id(**kwargs)

        # us that information to generate the name of the next motif atlas release id
        next_release = self.next_motif_release_id(**kwargs)

        print('Motif atlas parent release: %s, parent_index: %s, next_release: %s' % (parent,parent_index,next_release))
        self.logger.info('Motif atlas parent release: %s, parent_index: %s, next_release: %s' % (parent,parent_index,next_release))

        if parent == '0.0':
            parent = next_release

        current_nr_release_id = self.nr_release_id(**kwargs)

        nr_release_id_major = int(current_nr_release_id.split(".")[0])
        nr_release_id_minor = int(current_nr_release_id.split(".")[1])
        next_motif_atlas_release_major = int(next_release.split(".")[0])
        next_motif_atlas_release_minor = int(next_release.split(".")[1])

        print("Representative set release %s" % current_nr_release_id)
        self.logger.info("Representative set release %s" % current_nr_release_id)

        # only create a release when major numbers agree and representative
        # set minor number is four times the motif atlas minor number
        if nr_release_id_major == next_motif_atlas_release_major and next_motif_atlas_release_minor * 4 == nr_release_id_minor:
            return [Releases(
                current_nr_release_id,
                parent,
                parent_index,
                next_release,
            )]
        else:
            raise core.InvalidState("This is not a week to produce a motif atlas release")

    def has_data(self, *args, **kwargs):
        """
        Return False so the rest of this stage will run.
        """

        return False


    # def ifes_x_ray_only(self, nr_release_id, molecule_type):
    #     """
    #     For releases up to 3.99.

    #     Get a listing of all IFE's to use in clustering. The IFE's must be
    #     from the given list of structures, the ife must be representative for
    #     each class and the class should have the given resolution.
    #     The experimental method must be in MOTIF_ALLOWED_METHODS.
    #     This gives a subset of a representative set.

    #     molecule_type can be 'RNA' or 'DNA'

    #     :pdbs: The pdbs to get the best chains and models for.
    #     :returns: A dictionary mapping from pdb id to a set of the best chains
    #     and models.
    #     """

    #     if molecule_type == 'RNA':
    #         class_start = 'NR'
    #     elif molecule_type == 'DNA':
    #         class_start = 'DNA'

    #     with self.session() as session:
    #         # 2024-06-23 add filter on ifes.new_style
    #         classranks = mod.NrClassRank
    #         classes = mod.NrClasses
    #         ifes = mod.IfeInfo
    #         pdbs = mod.PdbInfo
    #         query = session.query(classranks).\
    #             join(classes, classes.name == classranks.nr_class_name).\
    #             join(ifes, ifes.ife_id == classranks.ife_id).\
    #             join(pdbs, pdbs.pdb_id == ifes.pdb_id).\
    #             filter(classranks.rank == 0).\
    #             filter(classes.nr_release_id == nr_release_id).\
    #             filter(classes.resolution == MOTIF_RESOLUTION_CUTOFF).\
    #             filter(classes.name.like('%s%%' % class_start)).\
    #             filter(ifes.new_style == True).\
    #             filter(pdbs.experimental_technique.in_(MOTIF_ALLOWED_METHODS)).\
    #             order_by(classranks.ife_id)

    #         if not query.count():
    #             raise core.InvalidState("No %s ifes found for nr %s" % (molecule_type,nr_release_id))

    #         ife_list = []
    #         for result in query:
    #             ife_id = result.ife_id
    #             if "6CZR" in ife_id:
    #                 # erroneously in a separate equivalence class
    #                 continue
    #             ife_list.append(ife_id)

    #         return ife_list


    def ifes(self, nr_release_id, molecule_type):
        """
        Allow EM structures with very good resolution.
        Allow X-ray structures with good resolution.

        Get a listing of all IFE's to use in clustering. The IFE's must be
        from the given list of structures, the ife must be representative for
        each class and the class should have the given resolution.
        The experimental method must be in MOTIF_ALLOWED_METHODS.
        This gives a subset of a representative set.

        molecule_type can be 'RNA' or 'DNA'

        :pdbs: The pdbs to get the best chains and models for.
        :returns: A dictionary mapping from pdb id to a set of the best chains
        and models.
        """

        # focus the search on ranks within these thresholds
        if molecule_type == 'RNA':
            em_class_start = 'NR_2.0'
            xray_class_start = 'NR_3.5'
        elif molecule_type == 'DNA':
            em_class_start = 'DNA_2.0'
            xray_class_start = 'DNA_3.5'

        # new for release 4.0: allow EM structures with res <= 2.0, cqs2 <= 9
        # new for release 4.0: restrict x-ray to cqs2 <= 15
        # new for release 4.0: take the highest-ranked EM or x-ray to meet the criteria
        with self.session() as session:
            classranks = mod.NrClassRank
            classes = mod.NrClasses
            cqs = mod.NrCqs
            ifes = mod.IfeInfo
            pdbs = mod.PdbInfo
            query = session.query(classranks).\
                join(classes, classes.name == classranks.nr_class_name).\
                join(ifes, ifes.ife_id == classranks.ife_id).\
                join(pdbs, pdbs.pdb_id == ifes.pdb_id).\
                join(cqs, cqs.ife_id == ifes.ife_id).\
                filter(classes.nr_release_id == nr_release_id).\
                filter(ifes.new_style == True).\
                filter(
                    or_(
                        and_(pdbs.experimental_technique == "ELECTRON MICROSCOPY", classes.name.like('%s%%' % em_class_start), cqs.cqs2 <= 9),
                        and_(pdbs.experimental_technique != "ELECTRON MICROSCOPY", classes.name.like('%s%%' % xray_class_start), cqs.cqs2 <= 15)
                    )
                ).\
                order_by(classranks.rank)

            if not query.count():
                raise core.InvalidState("No %s ifes found for nr %s" % (molecule_type,nr_release_id))

            ife_list = []
            handles = set()
            for result in query:
                # track the 5-digit code to identify an equivalence class
                handle = result.nr_class_name.split("_")[2].split(".")[0]
                if not handle in handles:
                    handles.add(handle)
                    ife_id = result.ife_id
                    if "6CZR" in ife_id:
                        # erroneously in a separate equivalence class
                        continue
                    ife_list.append(ife_id)
                    if not result.rank == 0:
                        self.logger.info('Using %s from %s with rank %d' % (ife_id,result.nr_class_name,result.rank))

            return ife_list


    def loops_to_exclude(self, **kwargs):
        """
        Load the loops to exclude. Sometimes for testing we want to
        exclude loops from analysis without updating the loop quality checks.
        This adds such functionality.
        """
        exclude = BLACKLIST
        exclude_file = kwargs.get('manual', {}).get('loop_exclude_file', None)
        if exclude_file:
            with open(exclude_file, 'rb') as raw:
                exclude.update(line.strip() for line in raw)
        return exclude


    def loop_position_border_units(self, loop_type, ifes, size_limit=None, **kwargs):
        """
        Get the list of loop ids to use in clustering, along with a
        mapping from position to border and unit id.
        These loops must be from IFE's in the given list and marked as
        valid in the loop quality step.
        This is where RNA loops with modified nucleotides *used to* get excluded.
        New on 2024-12-10:  allow status 1 or 3, to allow loops with
        modified nucleotides explicitly.

        Parameters
        ----------
        loop_type : str
            The type of loop to use, eg, IL, HL.
        ifes : list
            A list of ife ids to find loops in.

        Returns
        -------
        loops: str
            A list of loop ids to process.
        """

        exclude = self.loops_to_exclude(**kwargs)

        with self.session() as session:
            loops = mod.LoopInfo
            quality = mod.LoopQa
            pos = mod.LoopPositions
            ife_chains = mod.IfeChains
            chain_info = mod.ChainInfo
            units = mod.UnitInfo
            query = session.query(pos.loop_id,pos.position_2023,pos.position,pos.border,pos.unit_id,units.chain_index).\
                join(loops, loops.loop_id == pos.loop_id).\
                join(quality, quality.loop_id == loops.loop_id).\
                join(units, units.unit_id == pos.unit_id).\
                join(chain_info,
                     (chain_info.chain_name == units.chain) &
                     (chain_info.pdb_id == units.pdb_id)).\
                join(ife_chains,
                     ife_chains.chain_id == chain_info.chain_id).\
                filter(loops.deprecate == 0).\
                filter(quality.status.in_([1,3])).\
                filter(loops.type == loop_type).\
                filter(ife_chains.ife_id.in_(ifes)).\
                distinct()

            if size_limit is not None:
                self.logger.info("Applying loop size limit of %d" % size_limit)
                query = query.filter(loops.length < size_limit)

            # avoid IL with only 4 nucleotides; happens with strange chain breaks
            if loop_type == "IL":
                query = query.filter(loops.length > 4)

            # Execute the query and get the result
            results = query.all()

            print('Retrieved %s %s loop positions' % (len(results), loop_type))
            self.logger.info('Retrieved %s %s loop positions' % (len(results), loop_type))

            loop_position_to_border_unit_id = {}
            loop_id_to_border_count = {}
            # loop_ids_position_changed = set()
            for result in query:
                loop_id = result.loop_id
                if not loop_id in exclude:
                    if not loop_id in loop_position_to_border_unit_id:
                        loop_position_to_border_unit_id[loop_id] = {}
                        loop_id_to_border_count[loop_id] = 0
                    loop_position_to_border_unit_id[loop_id][result.position_2023] = (result.border,result.unit_id)
                    loop_id_to_border_count[loop_id] += result.border

                    # loop_id_to_model_symmetry_chain_index[loop_id]

                    # if not result.position == result.position_2023:
                    #     loop_position_to_border_unit_id[loop_id]['position'] = 'changed'
                    #     if not loop_id in loop_ids_position_changed:
                    #         self.logger.info("Position changed for loop %s" % loop_id)
                    #         loop_ids_position_changed.add(loop_id)

            # identify and remove loops with empty strands, for whatever reason
            if loop_type == 'HL':
                border_requirement = 2
            elif loop_type == 'IL':
                border_requirement = 4
            elif loop_type[0] == 'J':
                border_requirement = 2 * int(loop_type[1:])
            for loop_id in list(loop_position_to_border_unit_id.keys()):
                if loop_id_to_border_count[loop_id] < border_requirement:
                    self.logger.info("Excluding loop %s, only %d border positions" % (loop_id,loop_id_to_border_count[loop_id]))
                    for position, data in loop_position_to_border_unit_id[loop_id].items():
                        if not position == 'position':
                            self.logger.info("%s\t%2d\t%d\t%s" % (loop_id,position,data[0],data[1]))
                    del loop_position_to_border_unit_id[loop_id]

            # identify and remove loops generated by a single non-trivial symmetry operator
            symmetry_exclusion_counter = 0
            for loop_id in list(loop_position_to_border_unit_id.keys()):
                default_symmetry_found = False
                duplicate_chain_index = None
                symmetries = set()
                for position in loop_position_to_border_unit_id[loop_id].keys():
                    if not position == 'position':
                        unit_id = loop_position_to_border_unit_id[loop_id][position][1]
                        fields = unit_id.split("|")
                        if len(fields) < 9:
                            default_symmetry_found = True
                        else:
                            symmetries.add(fields[8])

                if not default_symmetry_found:
                    if len(symmetries) == 1:
                        symmetry = list(symmetries)[0]
                        if not symmetry == "P_1":
                            self.logger.info("Excluding loop %s, all symmetry operation %s" % (loop_id,symmetry))
                            del loop_position_to_border_unit_id[loop_id]
                            symmetry_exclusion_counter += 1

                if duplicate_chain_index:
                    self.logger.info("Excluding loop %s because of chain index duplication %s" % (loop_id,duplicate_chain_index))
                    del loop_position_to_border_unit_id[loop_id]

            self.logger.info("Excluding %d loops because of symmetries" % symmetry_exclusion_counter)

            # sort by length of first unit id to keep first strand without symmetries
            first_unit_id_length = []
            for loop_id in loop_position_to_border_unit_id.keys():
                if 1 in loop_position_to_border_unit_id[loop_id]:
                    unit_id = loop_position_to_border_unit_id[loop_id][1][1]
                    fields = unit_id.split("|")
                    first_unit_id_length.append((len(fields),loop_id))
                    self.logger.info('First unit id length for loop %s is %d' % (loop_id,len(fields)))

            # make sure we do not have the same loop twice, with different symmetries
            loop_short_unit_ids = set()
            for l, loop_id in sorted(first_unit_id_length):
                short_unit_id = []
                for position in loop_position_to_border_unit_id[loop_id].keys():
                    if not position == 'position':
                        unit_id = loop_position_to_border_unit_id[loop_id][position][1]
                        fields = unit_id.split("|")
                        short_unit_id.append("|".join(fields[:5]))
                short_unit_ids = ",".join(sorted(short_unit_id))
                if short_unit_ids in loop_short_unit_ids:
                    self.logger.info("Excluding loop %s because of duplicate short unit ids %s" % (loop_id,short_unit_ids))
                    del loop_position_to_border_unit_id[loop_id]
                else:
                    loop_short_unit_ids.add(short_unit_ids)

        return loop_position_to_border_unit_id

    def load_and_cache(self, loop_type, releases, folder):
        """
        This will load the result of running the pipeline and cache the
        data.

        :param str loop_type: The loop type, IL, HL, etc
        :param Release relases: The releases to use.
        :param str folder: The path to the results of the pipeline.
        :returns: Nothing
        """

        builder = Builder(self.config, self.session)
        motifs = builder(releases.parent, releases.current, folder)
        self.cache(loop_type, motifs)


    def cluster(self, loop_type, releases, ifes, molecule_type, **kwargs):
        """
        Cluster all motifs of the given loop type in the given ifes.

        :param str loop_type: The type of loop, ie IL or HL.
        :param Release releases: The named tuple of the type of loops.
        :param list ifes: The ifes to lookup loops in.
        :param kwargs: The keyword arguments. Notably the before_date key. This
        will change the date for the nr release and loop release used. It
        should be a date that a loop release and nr release exist for that
        date.
        :returns: None. All data is cached and nothing is returned.
        """

        cached_data = kwargs.get('manual', {}).get(loop_type, None)
        if cached_data:
            self.logger.info("Using cached data at %s" % cached_data)
            return self.cached(cached_data)

        size_limit = kwargs.get('manual', {}).get('loop_size_limit', None)
        if size_limit is not None:
            try:
                size_limit = int(size_limit)
            except:
                self.logger.error("Given size limit %s was not an int",
                                  size_limit)
                raise core.InvalidState("Bad size limit")

        #loops = self.loops(loop_type, ifes, size_limit=size_limit, **kwargs)

        loop_position_to_border_unit = self.loop_position_border_units(loop_type, ifes, size_limit=size_limit, **kwargs)

        if len(loop_position_to_border_unit) == 0:
            self.logger.info("No %s loops found, not clustering" % loop_type)
            return []

        self.logger.info("Found %d %s loops, starting to cluster them" % (len(loop_position_to_border_unit),loop_type))

        # The following two lines will call builder.py which uses ClusterMotif in cluster.py
        # which calls main in compare_and_cluster.py.
        builder = Builder(self.config, self.session)
        motifs = builder(loop_type, releases.parent, releases.current, loop_position_to_border_unit, molecule_type)

        self.cache(loop_type, motifs)
        self.logger.info("Done clustering %s" % loop_type)

        return motifs


    def add_group_id_to_diagnostic_file(self, loop_type, motifs):
        # change group number in a loop diagnostic into motif group id
        # this is not a particuarly good place to do this, but this is
        # one place where the variables are available

        # self.logger.info(motifs)

        cluster_number_to_group_id = {}
        for cluster_number, motif in enumerate(motifs['motifs']):
            cluster_number_to_group_id[cluster_number] = motif['motif_id']

        loop_annotation_file = os.path.join(motifs['motifs'][0]['2d'].split('2ds')[0], loop_type + '_similar_pairs.txt')

        self.logger.info('Modifying loop annotation file %s' % loop_annotation_file)

        if os.path.exists(loop_annotation_file):
            with open(loop_annotation_file, 'rt') as f:
                loop_annotations = f.readlines()

            # print('found:')
            # print(loop_annotations)

            with open(loop_annotation_file, 'wt') as f:
                for line in loop_annotations:
                    fields = line.strip().split('\t')
                    fields[2] = cluster_number_to_group_id.get(int(fields[2]), fields[2])
                    fields[8] = cluster_number_to_group_id.get(int(fields[8]), fields[8])
                    f.write('\t'.join(fields) + '\n')
                    print('\t'.join(fields))


    def data(self, release, **kwargs):
        """
        Compute the release for the loop types specified above.

        Later motif stages store the clustered data.

        :param list pdbs: The pdbs to store.
        :returns: A list of new releases, one for each type (eg IL, HL).
        """

        molecule_types = ['RNA']

        for molecule_type in molecule_types:
            # determine which ifes should be in the current release
            ifes = self.ifes(release.nr,molecule_type)
            self.logger.info("Found %d %s ifes in release %s" % (len(ifes),molecule_type,release.nr))

            # look up the date of the nr release this motif atlas release
            # is based on, to store in ml_releases under date
            # That is not the same as the date on which this code is running.
            nr_date = self.get_nr_release_date(release.nr)
            self.logger.info("Using date %s for release %s" % (nr_date,release.nr))

            # use loop types above, or specify manually
            if 'loop_type' in kwargs.get('manual', {}):
                loop_types = kwargs['manual']['loop_type'].split(',')
            else:
                loop_types = self.types

            # do all the loop types before returning any data or writing anything to the database
            data = []
            for loop_type in loop_types:
                # the cluster method stores the data in some cache,
                # which is then read in later stages and put in the database
                motifs = self.cluster(loop_type, release, ifes, molecule_type, **kwargs)

                # this registers that the clustering ran and the release was created
                # the data is not put into the database until all loop_types are done
                data.append(mod.MlReleases(ml_release_id=release.current,
                                        parent_ml_release_id=release.parent,
                                        date=nr_date,
                                        type=loop_type,
                                        index=release.parent_index + 1,
                                        nr_release_id=release.nr))

                if len(motifs) > 0:
                    self.add_group_id_to_diagnostic_file(loop_type,motifs)

        return data

