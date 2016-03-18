import os
import csv

from pymotifs import core
from pymotifs import utils
from pymotifs.utils import matlab
from pymotifs.models import LoopPositions
from pymotifs.models import LoopInfo
from pymotifs.units.info import Loader as UnitInfoLoader
from pymotifs.loops.extractor import Loader as InfoLoader
from pymotifs.loops.release import Loader as ReleaseLoader


class Loader(core.Loader):
    merge_data = True
    dependencies = set([UnitInfoLoader, InfoLoader, ReleaseLoader])
    allow_no_data = True

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)
        self.precomputed = self.config['locations']['loops_mat_files']

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(LoopInfo).filter_by(pdb_id=pdb)
            ids = [result.loop_id for result in query]

        if not ids:
            return True

        with self.session() as session:
            return session.query(LoopPositions).\
                filter(LoopPositions.loop_id.in_(ids)).\
                delete(synchronize_session=False)

        return True

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(LoopPositions).\
                join(LoopInfo, LoopInfo.loop_id == LoopPositions.loop_id).\
                filter(LoopInfo.pdb_id == pdb)
            return bool(query.count())

    def parse(self, filename):
        keys = ['loop_id', 'position', 'unit_id', 'bulge', 'flanking',
                'border']
        int_keys = ['position', 'bulge', 'flanking', 'border']
        data = []
        with open(filename, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for row in reader:
                entry = dict(zip(keys, row))
                for key in int_keys:
                    entry[key] = int(entry[key])
                data.append(entry)
            return data

    def known(self, pdb):
        mapping = {}
        with self.session() as session:
            query = session.query(LoopPositions).\
                join(LoopInfo, LoopInfo.loop_id == LoopPositions.loop_id).\
                filter(LoopInfo.pdb_id == pdb)

            for result in query:
                data = utils.row2dict(result)
                mapping[(result.loop_id, result.position)] = data

        return mapping

    def annotation_file(self, pdb):
        mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
        path = str(os.path.join(self.precomputed, pdb))
        try:
            if not os.path.exists(path):
                os.mkdir(path)
        except:
            raise core.InvalidState("Could not create %s for matlab" % path)

        [output_file, err_msg] = mlab.loadLoopPositions(path, nout=2)
        if err_msg != '':
            raise matlab.MatlabFailed(err_msg)

        return output_file

    def annotations(self, pdb, remove=True):
        output_file = self.annotation_file(pdb)
        data = self.parse(output_file)
        if remove:
            os.remove(output_file)
        return data

    def data(self, pdb, **kwargs):
        """Update loop_positions table by loading data from the mat files
        stored in the PrecomputedData folder
        """
        data = []
        known = self.known(pdb)
        for row in self.annotations(pdb):
            entry = known.get((row['loop_id'], row['position']), {})
            if not entry:
                self.logger.info("New loop %s found", row['loop_id'])

            entry.update(row)
            entry['loop_id'] = row['loop_id']
            print(entry)
            data.append(LoopPositions(**entry))

        return data
