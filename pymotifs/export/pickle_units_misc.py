"""
Export .pickle file of ions, ligands, and water unit centers
"""

import gzip
import numpy as np
import os
import pickle
from sqlalchemy import or_

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import DATA_FILE_DIRECTORY
from pymotifs.units.center_rotation import Loader as CenterRotationLoader


class Exporter(core.Loader):
    """
    """

    # General Setup
    compressed = False
    mark = False
    dependencies = set([CenterRotationLoader])


    def has_data(self, pdb, *args, **kwargs):
        # self.logger.info("has_data: pdb: %s" % str(pdb))
        filename = self.filename(pdb)
        # self.logger.info("has_data: filename: %s" % filename)
        if os.path.exists(filename) is True:
            self.logger.info("has_data: filename %s exists" % filename)
            return True
        self.logger.info("has_data: filename %s is missing" % filename)
        return False


    def remove():
        pass


    def filename(self,pdb,**kwargs):
        """
        Create the filename for the given PDB.

        Parameters
        ----------
        pdb :

        Returns
        -------
        filename : str
            The path to write to.
        """

        filename = pdb + '_misc.pickle.gz'

        return os.path.join(DATA_FILE_DIRECTORY,'units',filename)


    def to_process(self, pdbs, **kwargs):
        """
        Because it's hard to identify ligands, just loop over all given pdbs
        and if there are ions, ligands, or waters, write the .pickle file.
        """

        # get a list of all files in directory that end with _misc.pickle.gz
        directory = os.path.join(DATA_FILE_DIRECTORY,'units')
        files = [f for f in os.listdir(directory) if f.endswith('_misc.pickle.gz')]

        # extract the pdb ids from the files into a set
        written_pdb_ids = set([f.split('_')[0] for f in files])

        need_to_write = sorted(set(pdbs) - written_pdb_ids)

        return need_to_write


    def data(self, pdb, **kwargs):
        """
        Export .pickle file of ion, ligand, water unit centers
        """
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,
                            mod.UnitInfo.pdb_id,
                            mod.UnitInfo.unit_type_id,
                            mod.UnitCenters.name,
                            mod.UnitCenters.x,
                            mod.UnitCenters.y,
                            mod.UnitCenters.z).\
                            join(mod.UnitCenters, mod.UnitInfo.unit_id == mod.UnitCenters.unit_id).\
                            filter(mod.UnitInfo.pdb_id == pdb).\
                            filter(mod.UnitCenters.name == 'geometric')

            unit_id_to_data = {}
            for row in query:
                if not row.unit_id in unit_id_to_data:
                    unit_id_to_data[row.unit_id] = {}
                    unit_id_to_data[row.unit_id]['chain_index'] = 0
                unit_id_to_data[row.unit_id][row.name] = np.asarray([row.x,row.y,row.z])

        unit_ids = []
        centers = []

        for unit_id, data in unit_id_to_data.items():
            unit_ids.append(unit_id)
            centers.append(data['geometric'])

        return [unit_ids,centers]


    def process(self, pdb, **kwargs):
        """
        Write ion, ligand, water centers data for the given pdb file

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        filename = self.filename(pdb)

        pinfo = self.data(pdb)

        with gzip.open(filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

