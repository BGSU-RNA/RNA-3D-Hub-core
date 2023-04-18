"""Import all flanking interactions.

Runs fr3d on the given files to determine all flankng interactions and then
imports them into the database.
"""

import os
import csv


from pymotifs import core
from pymotifs.utils import matlab
from pymotifs import models as mod

from pymotifs.mat_files import Loader as MatLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.pdbs.info import Loader as PdbLoader


class Loader(core.SimpleLoader):
    """A loader to generate and import the interaction annotations for
    structures.
    """

    allow_no_data = True

    dependencies = set([MatLoader, UnitLoader, PdbLoader])

    @property
    def table(self):
        return mod.UnitPairsFlanking

    def query(self, session, pdb):
        """Create a query to access interaction data for the given pdb.

        :session: The database session to use.
        :pdb: The pdb id
        :returns: A query to get interaction data.
        """
        return session.query(mod.UnitPairsFlanking).filter_by(pdb_id=pdb)

    def parse(self, filename, pdb):
        """Reads the csv file, imports all interactions, deletes the file when
        done to avoid stale data and free up disk space

        :filename: The input filename.
        :pdb: The pdb id.
        :returns: A list of Interaction objects.
        """

        data = []
        with open(filename, 'rb') as raw:
            reader = csv.reader(raw, delimiter=',', quotechar='"')
            for index, row in enumerate(reader):
                if not row[0] or not row[1]:
                    msg = "Line %s of %s did not include both units"
                    raise core.InvalidState(msg % (index,filename))
                data.append({'unit_id_1':row[0], 'unit_id_2':row[1], 'flanking' : int(row[2]), 'pdb_id' : pdb})

        return data

    def data(self, pdb, **kwargs):
        """Compute the interaction annotations for a pdb file.

        :pdb: The pdb id to process.
        :kwargs: Keyword arguments.
        :returns: The interaction annotations.
        """
        mlab = matlab.Matlab(str(self.config['locations']['fr3d_root']))  # connect to Matlab
                                                                                #### https://github.com/BGSU-RNA/RNA-3D-Hub-core/blob/09d1044cd30bd396e701d6eb91b8eef75e78b1d4/conf/bootstrap.json.txt#L31
        self.logger.info('Running matlab on %s', pdb)
        ifn, status, err_msg = mlab.loadFlankings(pdb, nout=3)            # Matlab loads .mat file for this pdb and returns flanking pair list
        print(ifn)
        status = status[0][0]
        if status == 0:
            data = self.parse(ifn, pdb)
            os.remove(ifn)
            return data
        elif status == 2:
            raise core.Skip('PDB file %s has no nucleotides' % pdb)
        elif status == 3:
            raise core.Skip('PDB file %s has no flanking interactions' % pdb)
        raise core.InvalidState('Matlab error code %i when analyzing %s' %
                                status, pdb)
