"""
Import information about redundancy between nucleotides within a PDB file
"""

import os
import csv

from pymotifs import core
from pymotifs.models import UnitRedundancies


class RedundantNucleotidesLoader(core.SimpleLoader):

    def _parse(self, raw, pdb_file):
        data = []
        reader = csv.reader(raw, delimiter=',', quotechar='"')
        for row in reader:
            data.append({
                'unit_id1': row[0],
                'unit_id2': row[1],
                'pdb_id': pdb_file
            })
        return data

    def query(self, session, pdb):
        return session.query(UnitRedundancies).filter_by(pdb_id=pdb)

    def data(self, pdb):
        matlab = core.Matlab(self.config['locations']['fr3d_root'])
        ifn, err_msg = matlab.loadRedundantNucleotides(pdb, nout=2)
        if err_msg != '':
            raise core.MatlabFailed(err_msg)

        with open(ifn, 'rb') as raw:
            data = [UnitRedundancies(**unit) for unit in self._parse(raw, pdb)]
        os.remove(ifn)
        return data
