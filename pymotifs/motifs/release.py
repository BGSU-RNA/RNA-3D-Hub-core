"""This will create a new motif release and cluster all motifs. This will
cluster loops from all structure IFE's, that are representatives of their NR
class and come from X-ray structures. The clustered motifs are not yet stored
in the database but are just computed and temporary files are written.
"""

import datetime as dt

from sqlalchemy import desc

from pymotifs import core
from pymotifs.utils import releases as rel
from pymotifs import models as mod
from pymotifs.motifs.cluster import ClusterMotifs

from pymotifs.nr.release import Loader as NrReleaseLoader
from pymotifs.loops.release import Loader as LoopReleaseLoader

from pymotifs.nr.loader import Loader as NrLoader
from pymotifs.loops.loader import Loader as LoopLoader

""" A list of loop ids we should not use for clustering. """
BLACKLIST = set([
    'HL_3ICQ_004',
    'HL_3V2F_005',
    'HL_1Y0Q_002',
    'HL_2IL9_002',
    'HL_2IL9_005'
])


class Loader(core.MassLoader):
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

    def load_and_cache(self, folder):
        pass

    def data(self, pdbs, **kwargs):
        """Compute the releases for the given pdbs. This will cluster the
        motifs in the files as well as write a release. Future motif stages
        store the clustered data.

        :param list pdbs: The pdbs to store.
        :returns: A list of new releases, one for each type (eg IL, HL).
        """

        motif_release, index = self.current_id()
        loop_release = LoopReleaseLoader(self.config, self.session).\
            current_id()
        nr_release, _ = NrReleaseLoader(self.config, self.session).current_id()

        next = rel.next_id(motif_release,
                           mode=self.config['release_mode']['motifs'])

        self.logger.info("Current release is %s", motif_release)
        self.logger.info("Motif release will be %s", next)

        parent = motif_release
        if parent == '0.0':
            parent = next

        ifes = self.ifes(nr_release, pdbs)
        self.logger.info("Found %i ifes using release %s" %
                         (len(pdbs), nr_release))
        if not ifes:
            raise core.InvalidState("No ifes found for nr %s" % nr_release)

        cluster = ClusterMotifs(self.config, self.session.maker)

        now = dt.datetime.now()
        if kwargs.get('before_date', None):
            now = kwargs.get('before_date', dt.datetime.now())

        data = []
        for loop_type in self.types:
            folder = None
            if kwargs['dry_run']:
                self.logger.info("Skipping clustering in a dry run")
            else:
                loops = self.loops(loop_release, loop_type, ifes)
                self.logger.info("Found %i loops with release %s" %
                                 (len(loops), loop_release))

                if not loops:
                    raise core.InvalidState("No loops to cluster for %s" %
                                            loop_release)

                self.logger.info("Starting to cluster all %s", loop_type)
                folder = cluster(loop_type, loops)
                self.logger.info("Done clustering all %s into %s",
                                 loop_type, folder)
                self.load_and_cahe(folder)

            data.append(mod.MlReleases(id=next,
                                       parent_motif_release_id=parent,
                                       date=now,
                                       type=type,
                                       index=index + 1,
                                       loop_release_id=loop_release,
                                       nr_release_id=nr_release))

        return data
