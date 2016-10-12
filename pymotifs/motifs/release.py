"""This will create a new motif release and cluster all motifs. This will
cluster loops from all structure IFE's, that are representatives of their NR
class and come from X-ray structures. The clustered motifs are not yet stored
in the database but are just computed and temporary files are written.

Manual Options
--------------
nr_release : str
    The NR release to use, instead of the newest one.
loop_release : str
    The Loop release to use, instead of the newest one.
loop_size_limit : int
    The maximum loop size to use
loop_exclude_file : str
    Name of a file to that contains loops to exclude from clustering. This is
    meant for testing out different excluding procedures without update loop
    quality. The file should contain a single loop id per line, and only that.
"""

import datetime as dt
import collections as coll

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import releases as rel

from pymotifs.motifs.builder import Builder

from pymotifs.nr.release import Loader as NrReleaseLoader
from pymotifs.loops.release import Loader as LoopReleaseLoader

from pymotifs.nr.loader import Loader as NrLoader
from pymotifs.loops.loader import Loader as LoopLoader

Releases = coll.namedtuple('Releases', [
    'loop',
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
    resolution = '4.0'
    dependencies = set([NrLoader, LoopLoader])

    def nr_release_id(self, before_date=None, **kwargs):
        """Get the nr release. If no before_date is given then we get the
        latest, otherwise we get the release for the given date. If no release
        exists for that date then we fail.

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

    def loop_release_id(self, before_date=None, **kwargs):
        """Get the loop release. If no before_date is given then we get the
        latest, otherwise we get the release for the given date. If no release
        exists for that date then we fail.

        :param date before_date: The date to use.
        :returns: The loop release id.
        """

        if 'loop_release' in kwargs.get('manual', {}):
            return kwargs['manual']['loop_release']

        if before_date is None:
            return LoopReleaseLoader(self.config, self.session).current_id()

        with self.session() as session:
            query = session.query(mod.LoopReleases).\
                filter_by(date=before_date)

            if query.count() != 1:
                raise core.InvalidState("No loop release on %s", before_date)
            return query.one().loop_release_id

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

    def to_process(self, pdbs, **kwargs):
        """Determine what to should process. This stage works on a set of loop
        and nr releases.

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
            self.loop_release_id(**kwargs),
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
        each class and the class should have the given resolution. Also, the
        structure must be an x-ray structure.

        :pdbs: The pdbs to get the best chains and models for.
        :returns: A dictionary mapping from pdb id to a set of the best chains
        and models.
        """

        with self.session() as session:
            chains = mod.NrChains
            classes = mod.NrClasses
            ifes = mod.IfeInfo
            pdbs = mod.PdbInfo
            query = session.query(chains).\
                join(classes, classes.nr_class_id == chains.nr_class_id).\
                join(ifes, ifes.ife_id == chains.ife_id).\
                join(pdbs, pdbs.pdb_id == ifes.pdb_id).\
                filter(chains.rep == 1).\
                filter(chains.nr_release_id == nr_release_id).\
                filter(classes.resolution == self.resolution).\
                filter(pdbs.experimental_technique == 'X-RAY DIFFRACTION').\
                order_by(chains.ife_id)

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

    def loops(self, loop_release_id, loop_type, ifes, size_limit=None,
              **kwargs):
        """Get the list of loop ids to use in clustering. These loops must be
        from IFE's in the given list and marked as valid in the loop quality
        step.

        Parameters
        ----------
        loop_release_id : str
            The loop release id to use.
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
                filter(quality.loop_release_id == loop_release_id).\
                filter(loops.type == loop_type).\
                filter(ife_chains.ife_id.in_(ifes)).\
                filter(~loops.loop_id.in_(BLACKLIST)).\
                distinct()

            if size_limit is not None:
                query = query.filter(loops.length < size_limit)

            found.update(r.loop_id for r in query if r.loop_id not in exclude)

        if not loops:
            raise core.InvalidState("No loops to cluster for %s" %
                                    loop_release_id)
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

    def next_motif_release_id(self, **kwargs):
        """Get the next motif id.

        :returns: The next motif id.
        """

        motif_release, index = self.current_id()
        mode = self.config['release_mode']['motifs']
        return rel.next_id(motif_release, mode=mode)

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

        size_limit = kwargs.get('manual', {}).get('loop_size_limit', None)
        if size_limit is not None:
            try:
                size_limit = int(size_limit)
            except:
                self.logger.error("Give size limit %s was not an int",
                                  size_limit)
                raise core.InvalidState("Bad size limit")

        loops = self.loops(releases.loop, loop_type, ifes,
                           size_limit=size_limit, **kwargs)
        self.logger.info("Found %i loops", len(loops))
        self.logger.info("Starting to cluster all %s", loop_type)
        builder = Builder(self.config, self.session)
        motifs = builder(loop_type, releases.parent, releases.current, loops)
        self.cache(loop_type, motifs)
        self.logger.info("Done clustering %s", loop_type)

    def data(self, releases, **kwargs):
        """Compute the releases for the given pdbs. This will cluster the
        motifs in the files as well as write a release. Future motif stages
        store the clustered data.

        :param list pdbs: The pdbs to store.
        :returns: A list of new releases, one for each type (eg IL, HL).
        """

        ifes = self.ifes(releases.nr)
        self.logger.info("Found %i ifes", len(ifes))

        now = dt.datetime.now()
        if kwargs.get('before_date', None):
            now = kwargs.get('before_date', dt.datetime.now())

        data = []
        for loop_type in self.types:
            self.cluster(loop_type, releases, ifes, **kwargs)
            data.append(mod.MlReleases(ml_release_id=releases.current,
                                       parent_ml_release_id=releases.parent,
                                       date=now,
                                       type=loop_type,
                                       index=releases.parent_index + 1,
                                       loop_release_id=releases.loop,
                                       nr_release_id=releases.nr))

        return data
