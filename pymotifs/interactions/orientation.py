"""Import all anti/syn annotations.

Runs fr3d on the given files to determine nucleotide conformation and then
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
        return mod.UnitPairsOrientation

    def query(self, session, pdb):
        """Create a query to access interaction data for the given pdb.

        :session: The database session to use.
        :pdb: The pdb id
        :returns: A query to get interaction data.
        """
        return session.query(mod.UnitPairsOrientation).filter_by(pdb_id=pdb)

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
                '''
                The assumption here is that FR3D generates csv files with two columns:
                Column 1: unit_id
                Column 2: base_orientation (anti, syn or intermediate)
                '''
                data.append({'unit_id':row[0], 'orientation' : row[1], 'pdb_id' : pdb})

        return data

    def data(self, pdb, **kwargs):
        """Compute the anti/syn annotations for a pdb file.

        :pdb: The pdb id to process.
        :kwargs: Keyword arguments.
        :returns: The interaction annotations.
        """
        mlab = matlab.Matlab(str(self.config['locations']['fr3d_root']))

        self.logger.info('Running matlab on %s', pdb)
        # Need to change the matlab script name
        ifn, status, err_msg = mlab.loadFlankings(pdb, nout=3)
        status = status[0][0]
        if status == 0:
            data = self.parse(ifn, pdb)
            os.remove(ifn)
            return data
        # I'm not sure whether we need to change any of these    
        elif status == 2:
            raise core.Skip('PDB file %s has no nucleotides' % pdb)
        elif status == 3:
            raise core.Skip('PDB file %s has no flanking interactions' % pdb)
        raise core.InvalidState('Matlab error code %i when analyzing %s' %
                                status, pdb)
