from pymotifs.models import LoopInfo
from pymotifs import core
from pymotifs.utils import row2dict

from pymotifs.loops.extractor import Loader as LoopLoader


class Exporter(core.Exporter):
    headers = ['id', 'motif_id', 'pdb', 'nts']
    dependencies = set([LoopLoader])
    compressed = True

    def filename(self, pdb, **kwargs):
        return self.config['locations']['loops_gz']

    def loops(self, pdb):
        with self.session() as session:
            query = session.query(LoopInfo.loop_id.label('id'),
                                  LoopInfo.pdb_id.label('pdb'),
                                  LoopInfo.info.unit_ids.label('nts')
                                  ).\
                filter_by(pdb_id=pdb)

            count = query.count()
            if not count:
                self.logger.info("No loops found for %s", pdb)
            else:
                self.logger.info("Found %s loops for %s", count, pdb)

            return [row2dict(result) for result in query]

    def data(self, pdbs):
        for pdb in pdbs:
            self.logger.info("Writing out loops for %s", pdb)
            for loop in self.loops(pdb):
                yield loop
