from pymotifs import core

from pymotifs.models import MlLoops
from pymotifs.models import LoopInfo
from pymotifs.models import MlReleases

from pymotifs.utils import row2dict

from pymotifs.loops.extractor import Loader as LoopLoader
from pymotifs.loops.positions import Loader as PositionLoader


class Exporter(core.Exporter):
    headers = ['id', 'motif_id', 'pdb', 'nts']
    dependencies = set([LoopLoader, PositionLoader])
    compressed = True
    mark = False

    def filename(self, pdb, **kwargs):
        return self.config['locations']['loops_gz']

    def current_ml_release(self):
        with self.session() as session:
            current = session.query(MlReleases.ml_releases_id).\
                order_by(MlReleases.date).\
                limit(1).\
                first()

            if not current:
                return '0.0'

            return current.ml_release_id

    def loops(self, pdb):
        current_ml_release = self.current_ml_release()
        with self.session() as session:
            query = session.query(LoopInfo.loop_id.label('id'),
                                  LoopInfo.pdb_id.label('pdb'),
                                  LoopInfo.unit_ids.label('nts'),
                                  MlLoops.motif_id.label('motif_id')
                                  ).\
                outerjoin(MlLoops,
                          (MlLoops.loop_id == LoopInfo.loop_id) &
                          (MlLoops.ml_release_id == current_ml_release)).\
                filter(LoopInfo.pdb_id == pdb).\
                order_by(LoopInfo.loop_id)

            count = query.count()
            if not count:
                self.logger.info("No loops found for %s", pdb)
            else:
                self.logger.info("Found %s loops for %s", count, pdb)

            return [row2dict(result) for result in query]

    def data(self, pdbs, **kwargs):
        for pdb in pdbs:
            self.logger.info("Writing out loops for %s", pdb)
            for loop in self.loops(pdb):
                yield loop
