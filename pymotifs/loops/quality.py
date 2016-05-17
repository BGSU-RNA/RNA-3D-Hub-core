"""
Program for importing loop quality assurance data into the RNA 3D Hub database.

Some loops should be disqualified because they have unresolved or missing
nucleotides. See matlab code for disqualification codes.
"""

import os
import csv

from sqlalchemy import desc

from pymotifs import core
from pymotifs.utils import matlab
from pymotifs import models as mod

from pymotifs.loops.extractor import Loader as InfoLoader


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader])

    def current_id(self):
        with self.session() as session:
            query = session.query(mod.LoopReleases.loop_releases_id).\
                order_by(desc(mod.LoopReleases.date)).\
                limit(1)

            if query.count() == 0:
                return '0.0'

            return query.one().loop_releases_id

    def query(self, session, pdb):
        release_id = self.current_id()
        return session.query(mod.LoopQa).\
            join(mod.LoopInfo, mod.LoopInfo.loop_id == mod.LoopQa.loop_id).\
            filter(mod.LoopInfo.pdb_id == pdb).\
            filter(mod.LoopQa.release_id == release_id)

    def is_complementary(self, loop):
        pass

    def has_modified(self, loop):
        pass

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
                'loop_id': row[0],
                'status': row[1],
                'modifications': modres,
                'nt_signature': row[3],
                'complementary': compl,
                'release_id': release_id
            })

        return data

    def data(self, pdb, **kwargs):

        release_id = self.current_id()

        mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
        [ifn, err_msg] = mlab.aLoopQualityAssurance(pdb, nout=2)

        if err_msg != '':
            raise matlab.MatlabFailed('Error %s in pdb %s' % (err_msg, pdb))

        with open(ifn, 'rb') as raw:
            data = [mod.LoopQa(**loop) for loop in self.parse(raw, release_id)]

        os.remove(ifn)
        return data
