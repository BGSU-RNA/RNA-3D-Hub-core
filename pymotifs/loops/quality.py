"""
Program for importing loop quality assurance data into the RNA 3D Hub database.

Some loops should be disqualified because they have unresolved or missing
nucleotides. See matlab code for disqualification codes.
"""

import os
import csv

from pymotifs import core
from pymotifs.utils.releases import Release
from pymotifs.models import LoopQa
from pymotifs.models import LoopsAll
from pymotifs.loops.extractor import Loader as InfoLoader


class LoopQualityLoader(core.Loader):
    dependencies = set([InfoLoader])

    def has_data(self, pdb, **kwargs):
        helper = Release(self.config, self.session.maker)
        release_id = helper.current('loop')

        with self.session() as session:
            query = session.query(LoopQa.id).\
                join(LoopsAll, LoopsAll.id == LoopQa.id).\
                filter(LoopsAll.pdb == pdb).\
                filter(LoopQa.release_id == release_id)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        helper = Release(self.config, self.session.maker)
        release_id = helper.current('loop')
        with self.session() as session:
            query = session.query(LoopsAll.id).\
                filter_by(pdb=pdb, release_id=release_id)
            loop_ids = [result.id for result in query]

        if not loop_ids:
            return

        with self.session() as session:
            session.query(LoopQa).\
                filter(LoopQa.id.in_(loop_ids)).\
                delete(synchronize_session=False)

    def parse(self, raw, release_id):
        """Reads the csv file, imports all distances, deletes the file when done
           to avoid stale data and free up disk space
        """

        self.logger.info('Importing qa')
        reader = csv.reader(raw, delimiter=',', quotechar='"')
        data = []
        for i, row in enumerate(reader):
            modres = row[2]
            if modres == '':
                modres = None

            compl = row[4]
            if compl == '':
                compl = None

            data.append({
                'id': row[0],
                'status': row[1],
                'modifications': modres,
                'nt_signature': row[3],
                'complementary': compl,
                'release_id': release_id
            })

        return data

    def data(self, pdb, **kwargs):

        helper = Release(self.config, self.session.maker)
        release_id = helper.current('loop')

        mlab = core.Matlab(self.config['locations']['fr3d_root'])
        [ifn, err_msg] = mlab.aLoopQualityAssurance(pdb, nout=2)

        if err_msg != '':
            raise core.MatlabFailed('Error %s in pdb %s' % (err_msg, pdb))

        with open(ifn, 'rb') as raw:
            data = [LoopQa(**loop) for loop in self.parse(raw, release_id)]

        os.remove(ifn)
        return data
