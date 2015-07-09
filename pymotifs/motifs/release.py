import datetime

from fr3d import unit_ids as uid

from pymotifs import core
from pymotifs.motifs.cluster import ClusterMotifs
from pymotifs.utils.releases import Release
from pymotifs.models import PdbBestChainsAndModels
from pymotifs.models import LoopQa
from pymotifs.models import LoopsAll
from pymotifs.models import MlReleases
from pymotifs.models import NrPdbs
from pymotifs.models import NrClasses

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

    def best_chains(self, pdbs):
        """Get a mapping from pdb id to the best chains.

        :pdbs: The pdbs to get the best chains and models for.
        :returns: A dictonary mapping from pdb id to a set of the best chains
        and models.
        """

        best_chains = {}
        with self.session() as session:
            query = session.query(PdbBestChainsAndModels).\
                filter(PdbBestChainsAndModels.pdb_id.in_(pdbs))

            for result in query:
                best_chains[result.pdb_id] = set(result.best_chains.split(','))

        if not best_chains:
            raise core.StageFailed("Could not get best chains")

        return best_chains

    def valid_loops(self, loop_type, pdbs):
        """Get all valid loop ids that may go into clustering. A loop is valid
        if it is in the most recent release, is in the given list of pdbs, is
        of the requested type (IL, HL), have a valid QA status, and are not in
        BLACKLIST.
        """

        helper = Release(self.config, self.session.maker)
        release_id = helper.current('loop')

        with self.session() as session:
            query = session.query(LoopsAll.id, LoopsAll.nt_ids).\
                join(LoopQa, LoopQa.id == LoopsAll.id).\
                filter(LoopQa.status == 1).\
                filter(LoopQa.release_id == release_id).\
                filter(LoopsAll.type == loop_type).\
                filter(LoopsAll.pdb.in_(pdbs)).\
                filter(~(LoopsAll.id.in_(BLACKLIST)))

            data = []
            for result in query:
                data.append({
                    'id': result.id,
                    'unit_ids': result.nt_ids
                })
            return data

    def select_structures(self, pdbs):
        helper = Release(self.config, self.session.maker)
        nr_release = helper.current('nr')

        with self.session() as session:
            query = session.query(NrPdbs.id).\
                join(NrClasses, NrClasses.id == NrPdbs.class_id).\
                filter(NrClasses.release_id == NrPdbs.release_id).\
                filter(NrPdbs.release_id == nr_release).\
                filter(NrClasses.resolution == self.resolution).\
                filter(NrPdbs.rep == 1).\
                filter(NrPdbs.id.in_(pdbs))

            return [result.id for result in query]

    def loops(self, loop_type, pdbs):

        best_chains = self.best_chains(pdbs)
        valid = self.valid_loops(loop_type, pdbs)
        self.logger.info("Found %i valid loops", len(valid))

        """keep only loops from best chains based on their nt_ids"""
        loops = []
        for loop in valid:
            chains = set()
            pdb = None
            for unit_id in loop['unit_ids'].split(','):
                if '_' in unit_id:
                    parts = unit_id.split('_')
                    pdb = parts[0]
                    chains.add(parts[3])
                else:
                    data = uid.parse(unit_id)
                    pdb = data['pdb']
                    chains.add(data['chain'])

            if not chains:
                raise core.StageFailed("Could not determine loop chains")

            if chains.issubset(best_chains[pdb]):
                loops.append(loop['id'])

        self.logger.info('Selected %i loops', len(loops))

        return loops

    def data(self, pdbs, **kwargs):
        helper = Release(self.config, self.session.maker)
        current = helper.current('motif')
        next = helper.next(current, self.config['release_mode']['motifs'])
        self.logger.info("Motif Release will be %s", next)

        selected_pdbs = self.select_structures(pdbs)
        self.logger.info("Using %s of the %s structures",
                         len(selected_pdbs), pdbs)

        loops = self.loop(selected_pdbs)

        cluster = ClusterMotifs(self.config, self.session.maker)

        data = []
        for type in self.types:

            if kwargs['dry_run']:
                self.logger.info("Skipping clustering in a dry run")
                continue

            self.logger.info("Starting to cluster all %s", type)
            folder, annotation = cluster(type, loops)
            self.logger.info("Done clustering all %s", type)

            data.append(MlReleases(id=next, date=datetime.datetime.now(),
                                   type=type, description=folder,
                                   annotation=annotation))

        return data
