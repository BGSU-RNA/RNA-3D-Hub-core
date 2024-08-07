"""This will create a new motif release and cluster all motifs. This will
cluster loops from all structure IFE's, that are representatives of their NR
class and come from X-ray structures. The clustered motifs are not yet stored
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

import datetime as dt
import collections as coll

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import releases as rel

from pymotifs.constants import MOTIF_ALLOWED_METHODS
from pymotifs.constants import MOTIF_RESOLUTION_CUTOFF

from pymotifs.motifs.builder import Builder

from pymotifs.nr.release import Loader as NrReleaseLoader
from pymotifs.loops.release import Loader as LoopReleaseLoader

from pymotifs.nr.loader import Loader as NrLoader
from pymotifs.loops.loader import Loader as LoopLoader

Releases = coll.namedtuple('Releases', [
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
    'HL_2IL9_005'
])


class Loader(core.MassLoader):
    """A loader to create a motif release as well as process and cache all release
    data.
    """

    types = ['IL', 'HL']
    ## types = ['J3']
    dependencies = set([NrLoader, LoopLoader])

    def nr_release_id(self, before_date=None, **kwargs):
        """Get the nr release, if not given manually.
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
        """ Get the date of the specified nr release
        """

        with self.session() as session:
            query = session.query(mod.NrReleases).\
                filter_by(nr_release_id=nr_release_id)

            if query.count() != 1:
                raise core.InvalidState("Unclear date for nr release %s", nr_release_id)

            return query.one().date


    def current_id(self, **kwargs):
        """Get the current motif release id and the index. If there is no
        release then the release_id is 0.0 and the index is 0.

        :returns: The current release id and index.
        """

        with self.session() as session:
            query = session.query(mod.MlReleases.ml_release_id,
                                  mod.MlReleases.index).\
                order_by(desc(mod.MlReleases.index)).\
                limit(1)

            if query.count() == 0:
                return ('0.0', 0)

            current = query.one()
            return (current.ml_release_id, current.index)


    def next_motif_release_id(self, **kwargs):
        """Get the next motif id.

        :returns: The next motif id.
        """

        motif_release, index = self.current_id()
        mode = self.config['release_mode']['motifs']
        return rel.next_id(motif_release, mode=mode)


    def to_process(self, pdbs, **kwargs):
        """Determine what to process.
        This stage works on a set of nr releases.

        Returns
        -------
        releases : Release
            A Release object describing all releases to use.
        """

        parent, parent_index = self.current_id(**kwargs)
        next_release = self.next_motif_release_id(**kwargs)
        if parent == '0.0':
            parent = next_release

        return [Releases(
            self.nr_release_id(**kwargs),
            parent,
            parent_index,
            next_release,
        )]


    def has_data(self, *args, **kwargs):
        """This will always return True because we only want to update if the time
        difference has been large enough.
        """
        return False


    def ifes(self, nr_release_id):
        """Get a listing of all IFE's to use in clustering. The IFE's must be
        from the given list of structures, the ife must be representative for
        each class and the class should have the given resolution.
        The experimental method must be in MOTIF_ALLOWED_METHODS.
        This gives a subset of a representative set.

        :pdbs: The pdbs to get the best chains and models for.
        :returns: A dictionary mapping from pdb id to a set of the best chains
        and models.
        """

        # with self.session() as session:
        #     chains = mod.NrChains
        #     classes = mod.NrClasses
        #     ifes = mod.IfeInfo
        #     pdbs = mod.PdbInfo
        #     query = session.query(chains).\
        #         join(classes, classes.nr_class_id == chains.nr_class_id).\
        #         join(ifes, ifes.ife_id == chains.ife_id).\
        #         join(pdbs, pdbs.pdb_id == ifes.pdb_id).\
        #         filter(chains.rep == 1).\
        #         filter(chains.nr_release_id == nr_release_id).\
        #         filter(classes.resolution == MOTIF_RESOLUTION_CUTOFF).\
        #         filter(pdbs.experimental_technique.in_(MOTIF_ALLOWED_METHODS)).\
        #         order_by(chains.ife_id)

        ## made a new query because we stop updating the nr_chains table
        with self.session() as session:
            query = session.query(mod.NrClassRank.ife_id).\
                join(mod.NrClasses, mod.NrClasses.name == mod.NrClassRank.nr_class_name).\
                join(mod.IfeInfo, mod.IfeInfo.ife_id == mod.NrClassRank.ife_id).\
                join(mod.PdbInfo,mod.PdbInfo.pdb_id == mod.IfeInfo.pdb_id).\
                filter(mod.NrClassRank.rank == 0).\
                filter(mod.NrClasses.nr_release_id == nr_release_id).\
                filter(mod.NrClasses.resolution == MOTIF_RESOLUTION_CUTOFF).\
                filter(mod.PdbInfo.experimental_technique.in_(MOTIF_ALLOWED_METHODS)).\
                order_by(mod.NrClassRank.ife_id)

                # omitted this for now to run the motif atlas release with old ifes
                # filter(mod.IfeInfo.new_style == True).\



            if not query.count():
                raise core.InvalidState("No ifes found for nr %s" %
                                        nr_release_id)

            return [result.ife_id for result in query]


    def loops_to_exclude(self, **kwargs):
        """Load the loops to exclude. Sometimes for testing usages we want to
        exclude loops from analysis without updating the loop quality checks.
        This adds such functionality.
        """
        exclude = BLACKLIST
        exclude_file = kwargs.get('manual', {}).get('loop_exclude_file', None)
        if exclude_file:
            with open(exclude_file, 'rb') as raw:
                exclude.update(line.strip() for line in raw)
        return exclude


    def loops(self, loop_type, ifes, size_limit=None,
              **kwargs):
        """
        Deprecated code as of 2022-10-23 since it cannot deal with symmetry operators.

        Get the list of loop ids to use in clustering. These loops must be
        from IFE's in the given list and marked as valid in the loop quality
        step.
        This is where RNA loops with modified nucleotides get excluded.

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

        found = set()
        with self.session() as session:
            loops = mod.LoopInfo
            quality = mod.LoopQa
            pos = mod.LoopPositions
            ife_chains = mod.IfeChains
            chain_info = mod.ChainInfo
            units = mod.UnitInfo
            query = session.query(loops.loop_id).\
                join(quality, quality.loop_id == loops.loop_id).\
                join(pos, pos.loop_id == loops.loop_id).\
                join(units, units.unit_id == pos.unit_id).\
                join(chain_info,
                     (chain_info.chain_name == units.chain) &
                     (chain_info.pdb_id == units.pdb_id)).\
                join(ife_chains,
                     ife_chains.chain_id == chain_info.chain_id).\
                filter(quality.status == 1).\
                filter(loops.type == loop_type).\
                filter(ife_chains.ife_id.in_(ifes)).\
                filter(~loops.loop_id.in_(BLACKLIST)).\
                distinct()

            if size_limit is not None:
                query = query.filter(loops.length < size_limit)

            # avoid IL with only 4 nucleotides; happens with strange chain breaks
            if loop_type == "IL":
                query = query.filter(loops.length > 4)


            found.update(r.loop_id for r in query if r.loop_id not in exclude)

        if not found:
            raise core.InvalidState("No loops to cluster for %s" %
                                    loop_release_id)
        return sorted(found)


    def loops_and_strands(self, loop_type, ifes, size_limit=None,
              **kwargs):
        """
        Get the list of loop ids to use in clustering, along with a
        list of unit ids for each strand.
        These loops must be from IFE's in the given list and marked as
        valid in the loop quality step.
        This is where RNA loops with modified nucleotides get excluded.

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

        found = set()
        with self.session() as session:
            loops = mod.LoopInfo
            quality = mod.LoopQa
            pos = mod.LoopPositions
            ife_chains = mod.IfeChains
            chain_info = mod.ChainInfo
            units = mod.UnitInfo
            query = session.query(pos.loop_id,pos.position,pos.border,pos.unit_id).\
                join(loops, loops.loop_id == pos.loop_id).\
                join(quality, quality.loop_id == loops.loop_id).\
                join(units, units.unit_id == pos.unit_id).\
                join(chain_info,
                     (chain_info.chain_name == units.chain) &
                     (chain_info.pdb_id == units.pdb_id)).\
                join(ife_chains,
                     ife_chains.chain_id == chain_info.chain_id).\
                filter(quality.status == 1).\
                filter(loops.type == loop_type).\
                filter(ife_chains.ife_id.in_(ifes)).\
                filter(~pos.loop_id.in_(BLACKLIST)).\
                distinct()

            if size_limit is not None:
                query = query.filter(loops.length < size_limit)

            # avoid IL with only 4 nucleotides; happens with strange chain breaks
            if loop_type == "IL":
                query = query.filter(loops.length > 4)

            loopdata = {}
            for result in query:
                loop_id = result.loop_id
                if not loop_id in exclude:
                    if not loop_id in loopdata:
                        loopdata[loop_id] = {}
                    loopdata[loop_id][result.position] = (result.border,result.unit_id)

            # identify and remove loops generated by a single non-trivial symmetry operator
            symmetry_exclusion_counter = 0
            for loop_id in list(loopdata.keys()):
                default_symmetry_found = False
                symmetries = set()
                for position in loopdata[loop_id].keys():
                    unit_id = loopdata[loop_id][position][1]
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
                            del loopdata[loop_id]
                            symmetry_exclusion_counter += 1

            self.logger.info("Excluding %d loops because of symmetries" % symmetry_exclusion_counter)

            self.logger.info(loopdata)

            found = loopdata.keys()

        if not found:
            raise core.InvalidState("No loops to cluster for %s" % loop_release_id)

        return sorted(found)

    def load_and_cache(self, loop_type, releases, folder):
        """This will load the result of running the pipeline and cache the
        data.

        :param str loop_type: The loop type, IL, HL, etc
        :param Release relases: The releases to use.
        :param str folder: The path to the results of the pipeline.
        :returns: Nothing
        """

        builder = Builder(self.config, self.session)
        motifs = builder(releases.parent, releases.current, folder)
        self.cache(loop_type, motifs)


    def cluster(self, loop_type, releases, ifes, **kwargs):
        """This will cluster all motifs of the given loop type in the given
        ifes. It will run the matlab code to cluster all motifs and then the
        produced data will be loaded into a single data structure and cached
        for future use.

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

        loops = self.loops_and_strands(loop_type, ifes, size_limit=size_limit, **kwargs)

        self.logger.info("Found %d %s loops" % (len(loops),loop_type))
        self.logger.info("Starting to cluster all %s" % loop_type)
        builder = Builder(self.config, self.session)
        motifs = builder(loop_type, releases.parent, releases.current, loops)
        self.cache(loop_type, motifs)
        self.logger.info("Done clustering %s" % loop_type)
        return motifs

    def data(self, releases, **kwargs):
        """Compute the releases for the given pdbs. This will cluster the
        motifs in the files as well as write a release. Future motif stages
        store the clustered data.

        :param list pdbs: The pdbs to store.
        :returns: A list of new releases, one for each type (eg IL, HL).
        """

        ifes = self.ifes(releases.nr)
        nr_date = self.get_nr_release_date(releases.nr)

        self.logger.info("Found %d ifes in release %s" % (len(ifes),releases.nr))
        self.logger.info("Using date %s for release %s" % (nr_date,releases.nr))


        """
        if ":" in now:
            nowstring = dt.strftime(dt.strptime(now, "%Y-%m-%d %H:%M:%S").date(), "%Y%m%d")
        else:
            nowstring = now.replace("-","")
        """

        loop_types = self.types

        if 'loop_type' in kwargs.get('manual', {}):
            loop_types = [kwargs['manual']['loop_type']]

        data = []
        for loop_type in loop_types:
            motifs = self.cluster(loop_type, releases, ifes, **kwargs)
            with open(motifs['graph'], 'rb') as raw:
                graphml = raw.read().replace('\n', '')
            data.append(mod.MlReleases(ml_release_id=releases.current,
                                       parent_ml_release_id=releases.parent,
                                       date=nr_date,
                                       type=loop_type,
                                       index=releases.parent_index + 1,
                                       nr_release_id=releases.nr,
                                       graphml=graphml))

        return data
