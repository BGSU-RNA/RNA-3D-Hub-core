"""This will create a new motif release and cluster all motifs. This will
cluster loops from all structure IFE's, that are representatives of their NR
class and come from X-ray structures. The clustered motifs are not yet stored
in the database but are just computed and temporary files are written.
"""

import datetime as dt
import collections as coll

from sqlalchemy import desc

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import releases as rel

from pymotifs.motifs.builder import Builder
from pymotifs.motifs.builder import MutualDiscrepancyLoader

from pymotifs.motifs.cluster import ClusterMotifs

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

"""The filename to cache all loop discrepancy data in"""
DISC_CACHE = 'disc'


class Loader(core.MassLoader):
    """A loader to create a motif release as well as process and cache all release
    data.
    """

    types = ['IL', 'HL']
    resolution = '4.0'
    dependencies = set([NrLoader, LoopLoader])

    def has_data(self, *args, **kwargs):
        """This will always return True because we only want to update if the time
        difference has been large enough.
        """
        return False

    def ifes(self, nr_release_id, pdb_ids):
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
                filter(ifes.pdb_id.in_(pdb_ids)).\
                filter(ifes.has_structured == 1).\
                filter(pdbs.experimental_technique == 'X-RAY DIFFRACTION').\
                order_by(chains.ife_id)

            if not query.count():
                raise core.InvalidState("No ifes found for nr %s" %
                                        nr_release_id)

            return [result.ife_id for result in query]

    def loops(self, loop_release_id, loop_type, ifes):
        """Get the list of loop ids to use in clustering. These loops must be
        from IFE's in the given list and marked as valid in the loop quality
        step.

        :param str loop_release_id: The loop release id to use.
        :param str loop_type: The type of loop to use, eg, IL, HL.
        :param list ifes: A list of ifes to find loops in.
        """

        found = []
        for ife in ifes:
            with self.session() as session:
                loops = mod.LoopInfo
                quality = mod.LoopQa
                pos = mod.LoopPositions
                ifes = mod.IfeInfo
                ife_chains = mod.IfeChains
                chain_info = mod.ChainInfo
                units = mod.UnitInfo
                query = session.query(loops.loop_id).\
                    join(quality, quality.loop_id == loops.loop_id).\
                    join(pos, pos.loop_id == loops.loop_id).\
                    join(units, units.unit_id == pos.unit_id).\
                    join(chain_info, chain_info.chain_name == units.chain).\
                    join(ife_chains,
                         ife_chains.chain_id == chain_info.chain_id).\
                    filter(quality.status == 1).\
                    filter(quality.loop_release_id == loop_release_id).\
                    filter(loops.type == loop_type).\
                    filter(ife_chains.ife_id == ife).\
                    distinct()

                found.extend(result.loop_id for result in query)

        if not loops:
            raise core.InvalidState("No loops to cluster for %s" %
                                    loop_release_id)
        return found

    def current_id(self):
        """Get the current motif release id and the index. If there is no
        release then the release_id is 0.0 and the index is 0.

        :returns: The current release id and index.
        """

        with self.session() as session:
            query = session.query(mod.MlReleases.ml_releases_id,
                                  mod.MlReleases.index).\
                order_by(desc(mod.MlReleases.index)).\
                limit(1)

            if query.count() == 0:
                return '0.0', 0

            current = query.one()
            return current.nr_release_id, current.index

    def load_and_cache(self, type, releases, folder):
        """This will load the result of running the pipeline and cache the
        data.

        :param str type: The loop type, IL, HL, etc
        :param Release relases: The releases to use.
        :param str folder: The path to the results of the pipeline.
        :returns: Nothing
        """

        builder = Builder(self.config, self.session)
        motifs = builder(releases.parent, releases.current, folder)
        self.cache(motifs, type)

    def cache_discrepancies(self, loop_type, folder):
        """This will load and cache discrepancies. If we have already cached
        some data we will append to the cached data as well.

        :param str loop_type: The loop type, IL, HL, etc.
        :param str folder: The path to where the motif data was computed.
        """

        loader = MutualDiscrepancyLoader(self.config, self.session)
        data = loader(folder)
        current = self.cached('disc')
        if current is not None:
            current.append(data)
        self.cache(data, DISC_CACHE)

    def nr_release_id(self, before_date=None, **kwargs):
        """Get the nr release. If no before_date is given then we get the
        latest, otherwise we get the release for the given date. If no release
        exists for that date then we fail.

        :param date before_date: The date to use.
        :returns: The nr release id.
        """

        if before_date is None:
            nr_release, _ = NrReleaseLoader(self.config, self.session).\
                current_id()
            return nr_release

        with self.session() as session:
            query = session.query(mod.NrReleases).\
                filter_by(date=before_date)

            if query.count() != 1:
                raise core.InvalidState("No nr release for date %s",
                                        before_date)
            return query.one().nr_release_id

    def loop_release_id(self, before_date=None, **kwargs):
        """Get the loop release. If no before_date is given then we get the
        latest, otherwise we get the release for the given date. If no release
        exists for that date then we fail.

        :param date before_date: The date to use.
        :returns: The loop release id.
        """

        if before_date is None:
            return LoopReleaseLoader(self.config, self.session).current_id()

        with self.session() as session:
            query = session.query(mod.LoopReleases).\
                filter_by(date=before_date)

            if query.count() != 1:
                raise core.InvalidState("No loop release for date %s",
                                        before_date)
            return query.one().loop_release_id

    def motif_relase_id(self, **kwargs):
        """Get the next motif id.

        :returns: The next motif id.
        """

        motif_release, index = self.current_id()
        parent = motif_release
        mode = self.config['release_mode']['motifs']
        current = rel.next_id(motif_release, mode=mode)

        if parent == '0.0':
            parent = current

        return current

    def releases(self, **kwargs):
        """Compute all releases.

        :returns: A Release object describing all releases to use.
        """

        return Releases(
            self.loop_release_id(**kwargs),
            self.nr_release_id(**kwargs),
            *self.loop_release_id(**kwargs)
        )

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

        cluster = ClusterMotifs(self.config, self.session.maker)
        loops = self.loops(releases.loop, loop_type, ifes)
        self.logger.info("Found %i loops", len(loops))

        self.logger.info("Starting to cluster all %s", loop_type)
        folder = cluster(loop_type, loops)
        self.logger.info("Done clustering all %s into %s",
                         loop_type, folder)
        self.load_and_cache(loop_type, releases, folder)
        self.cache_discrepancies(loop_type, folder)

    def data(self, pdbs, **kwargs):
        """Compute the releases for the given pdbs. This will cluster the
        motifs in the files as well as write a release. Future motif stages
        store the clustered data.

        :param list pdbs: The pdbs to store.
        :returns: A list of new releases, one for each type (eg IL, HL).
        """

        releases = self.releases(**kwargs)
        self.logger.info("Current releases are %s", str(releases))

        ifes = self.ifes(releases.nr, pdbs)
        self.logger.info("Found %i ifes" % (len(ifes)))

        now = dt.datetime.now()
        if kwargs.get('before_date', None):
            now = kwargs.get('before_date', dt.datetime.now())

        data = []
        for loop_type in self.types:
            self.cluster(loop_type, releases, ifes, **kwargs)
            data.append(mod.MlReleases(id=releases.current,
                                       parent_motif_release_id=releases.parent,
                                       date=now,
                                       type=loop_type,
                                       index=releases.parent_index + 1,
                                       loop_release_id=releases.loop,
                                       nr_release_id=releases.nr))

        return data
