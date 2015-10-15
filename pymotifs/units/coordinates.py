import os
import csv

import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitCoordinates
from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader])

    def query(self, session, pdb):
        return session.query(UnitCoordinates).\
            join(UnitInfo, UnitInfo.unit_id == UnitCoordinates.unit_id).\
            filter(UnitInfo.pdb_id == pdb)

    def parse(self, filename):
        data = []
        with open(filename, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for row in reader:
                data.append(UnitCoordinates(id=row[0], coordinates=row[9]))
        return data

    def data(self, pdb, **kwargs):
        mlab = core.Matlab()
        ifn, status, err_msg = mlab.loadCoordinates(pdb, nout=3)
        status = status[0][0]
        if status == 0:
            data = self.parse(ifn)
            os.remove(ifn)
            return data
        elif status == 2:
            self.logger.warning('File %s has no nucleotides', pdb)
        else:
            raise core.MatlabFailed('Matlab error code %i when analyzing %s' %
                                    (status, pdb))
