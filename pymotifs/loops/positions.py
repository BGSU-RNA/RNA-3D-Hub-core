import os
import csv

from pymotifs import core
from pymotifs import utils
from pymotifs.models import LoopPositions
from pymotifs.models import LoopsAll


class Loader(core.SimpleLoader):
    merge_data = True

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)
        self.matlab = core.Matlab(self.config['locations']['fr3d_root'])
        self.precomputed = self.config['locations']['loops_mat_files']

    def parse(self, filename):
        with open(filename, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for row in reader:
                yield row

    def has_data(self, pdb, **kwargs):
        return False

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            session.query(LoopPositions).\
                join(LoopsAll, LoopsAll.id == LoopPositions.loop_id).\
                filter(LoopsAll.pdb == pdb).\
                delete(synchronize_session='fetch')

    def known(self, pdb):
        mapping = {}
        with self.session() as session:
            query = session.query(LoopPositions).\
                join(LoopsAll, LoopsAll.id == LoopPositions.loop_id).\
                filter(LoopsAll.pdb == pdb)

            for result in query:
                data = utils.row2dict(result)
                mapping[(result.loop_id, result.position)] = data

        return mapping

    def get_annotations(self, pdb):
        path = os.path.join(self.precomputedData, pdb)
        [output_file, err_msg] = self.mlab.loadLoopPositions(path, nout=2)
        if err_msg != '':
            raise core.MatlabFailed(err_msg)

        return output_file

    def data(self, pdb):
        """Update loop_positions table by loading data from the mat files
        stored in the PrecomputedData folder
        """
        data = []
        output_file = self.get_annotations(pdb)
        known = self.known(pdb)
        for row in self.parse(output_file):
            (loop_id, position, nt_id, bulge, flanking, border) = row
            entry = known.get((loop_id, position), {})
            if not entry:
                self.logger.info("New loop %s found", loop_id)

            entry.update({
                'loop_id': loop_id,
                'position': position,
                'nt_id': nt_id,
                'flanking': int(flanking),
                'bulge': int(bulge),
                'border': int(border)
            })
            data.append(LoopPositions(**entry))

        os.remove(output_file)  # delete temporary csv file

        return data
