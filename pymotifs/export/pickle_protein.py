"""
Export .pickle file of protein unit centers
"""

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
        self.logger.info("has_data: pdb: %s" % str(pdb))
        filename = self.filename(pdb)
        self.logger.info("has_data: filename: %s" % filename)
        if os.path.exists(filename) is True:
            self.logger.info("has_data: filename %s exists" % filename)
            return True
        self.logger.info("has_data: filename %s is missing" % filename)
        return False


    def remove():
        pass


    def filename(self,pdb,**kwargs):
        """Create the filename for the given PDB.

        Parameters
        ----------
        pdb :

        Returns
        -------
        filename : str
            The path to write to.
        """

        filename = pdb + '_protein.pickle'

        return os.path.join(DATA_FILE_DIRECTORY,filename)


    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitInfo.pdb_id).\
                            filter(mod.UnitInfo.unit_type_id == 'aa').distinct()
        pdb_ids = []
        for row in query:
            pdb_ids.append(row.pdb_id)
        return pdb_ids


    def data(self, pdb, **kwargs):
        """
        Export .pickle file of protein unit centers
        """
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,
                            mod.UnitInfo.pdb_id,
                            mod.UnitInfo.chain_index,
                            mod.UnitInfo.unit_type_id,
                            mod.UnitCenters.name,
                            mod.UnitCenters.x,
                            mod.UnitCenters.y,
                            mod.UnitCenters.z).\
                            join(mod.UnitCenters, mod.UnitInfo.unit_id == mod.UnitCenters.unit_id).\
                            filter(mod.UnitInfo.pdb_id==pdb).\
                            filter(or_(mod.UnitCenters.name == 'aa_fg',mod.UnitCenters.name == 'aa_backbone')).\
                            order_by(mod.UnitInfo.chain_index)

            unit_id_to_data = {}
            for row in query:
                if not row.unit_id in unit_id_to_data:
                    unit_id_to_data[row.unit_id] = {}
                    unit_id_to_data[row.unit_id]['chain_index'] = row.chain_index
                if row.name == 'aa_fg':
                    unit_id_to_data[row.unit_id]['aa_fg'] = np.asarray([row.x,row.y,row.z])
                elif row.name == 'aa_backbone':
                    unit_id_to_data[row.unit_id]['aa_backbone'] = np.asarray([row.x,row.y,row.z])

        unit_ids = []
        chain_indices = []
        aa_fg = []
        aa_backbone = []

        for unit_id, data in unit_id_to_data.items():
            unit_ids.append(unit_id)
            chain_indices.append(data['chain_index'])
            if 'aa_fg' in data:
                aa_fg.append(data['aa_fg'])
            else:
                aa_fg.append(np.asarray([np.NaN,np.NaN,np.NaN]))
            if 'aa_backbone' in data:
                aa_backbone.append(data['aa_backbone'])
            else:
                aa_backbone.append(np.asarray([np.NaN,np.NaN,np.NaN]))

        # old method, which was fragile
        # for row in query:
        #     if row.unit_id not in ids:
        #         ids.append(row.unit_id)
        #         chain_indices.append(int(str(row.chain_index).replace('L','')))
        #         if row.name == 'aa_fg':
        #             aa_fg.append((np.asarray([row.x,row.y,row.z])))
        #             aa_backbone.append(np.asarray([np.NaN,np.NaN,np.NaN]))
        #         elif row.name == 'aa_backbone':
        #             aa_backbone.append((np.asarray([row.x,row.y,row.z])))
        #             aa_fg.append(np.asarray([np.NaN,np.NaN,np.NaN]))
        #     else:
        #         if row.name == 'aa_fg':
        #             aa_fg[len(aa_fg)-1]=(np.asarray([row.x,row.y,row.z]))
        #         elif row.name == 'aa_backbone':
        #             aa_backbone[len(aa_fg)-1]=(np.asarray([row.x,row.y,row.z]))


        return [unit_ids,chain_indices,aa_fg,aa_backbone]


    def process(self, pdb, **kwargs):
        """
        Write amino acid centers data for the given pdb file

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        filename = self.filename(pdb)

        pinfo = self.data(pdb)

        with open(filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

