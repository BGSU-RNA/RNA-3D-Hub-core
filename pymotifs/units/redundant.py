"""
Import information about redundancy between nucleotides within a PDB file
"""

import os
import csv

from pymotifs import core
from pymotifs.utils import matlab
from pymotifs import models as mod
from pymotifs.mat_files import Loader as MatLoader


class RedundantNucleotidesLoader(core.SimpleLoader):
    dependencies = set([MatLoader])
    allow_no_data = True
    table = mod.UnitRedundancies

    def _parse(self, raw, pdb_file):
        data = []
        reader = csv.reader(raw, delimiter=',', quotechar='"')
        for row in reader:
            data.append({
                'unit_id_1': row[0],
                'unit_id_2': row[1],
                'pdb_id': pdb_file
            })
        return data

    def query(self, session, pdb):
        query = session.query(mod.UnitRedundancies).filter_by(pdb_id=pdb)
        return query

    def data(self, pdb, **kwargs):
        mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
        ifn, err_msg = mlab.loadRedundantNucleotides(pdb, nout=2)
        if err_msg != '':
            raise matlab.MatlabFailed(err_msg)

        with open(ifn, 'rb') as raw:
            data = self._parse(raw, pdb)
        os.remove(ifn)
        return data
