import csv
import gzip
import shutil
import tempfile

from pymotifs import core
from pymotifs.models import UnitPairsInteractions
from pymotifs.interactions.pairwise import Loader as InterLoader


class InteractionExporter(core.Stage):
    headers = ['unit_id1', 'unit_id2', 'FR3D basepair (f_lwbp)',
               'FR3D stacking (f_stacks)', 'FR3D base phosphate (f_bphs)']
    dependencies = set([InterLoader])

    def is_missing(self, *args, **kwargs):
        return True

    def remove(self, *args, **kwargs):
        pass

    def compress(self, temp, **kwargs):
        temp_output_file = self.filename() + '-temp'
        temp.seek(0)
        handle = gzip.open(temp_output_file, 'wb')
        handle.writelines(temp)
        handle.close()
        shutil.move(temp_output_file, self.filename())

    def interactions(self, pdb):
        data = []
        with self.session() as session:
            query = session.query(UnitPairsInteractions).filter_by(pdb_id=pdb)
            for result in query:
                data.append({
                    'unit_id1': result.unit1_id,
                    'unit_id2': result.unit2_id,
                    'FR3D basepair (f_lwbp)': result.f_lwbp,
                    'FR3D stacking (f_stacks)': result.f_stacks,
                    'FR3D base phosphate (f_bphs)': result.f_bphs
                })
        self.logger.debug("Found %s interactions", len(data))
        return data

    def filename(self, *args, **kwargs):
        return self.config['locations']['interactions_gz']

    def process(self, pdb_id, **kwargs):
        self.logger.info('Writing out interactions for %s' % pdb_id)
        interactions = self.interactions(pdb_id)
        self.writer.writerows(interactions)

    def __call__(self, *args, **kwargs):
        with tempfile.TemporaryFile() as temp:
            self.writer = csv.DictWriter(temp, self.headers, quotechar='"',
                                         quoting=csv.QUOTE_ALL)
            self.writer.writerow(dict(zip(self.headers, self.headers)))

            super(InteractionExporter, self).__call__(*args, **kwargs)

            self.compress(temp, **kwargs)
