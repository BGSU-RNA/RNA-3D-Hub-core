"""Module for export of resolution and methods of existing pdb files
"""

import numpy as np
import os
import pickle

from pymotifs import core
from pymotifs import models as mod

from pymotifs.pdbs.info import Loader as InfoLoader

from collections import defaultdict

from os import path

from sqlalchemy import or_
from sqlalchemy.orm import aliased
from sqlalchemy.sql import select
from sqlalchemy.sql import union



class Exporter(core.Loader):
    """Export pairs data in pickle format.
    """


    # General Setup
    compressed = False 
    mark = False 
    dependencies = set([InfoLoader])
    #dependencies = set()



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


        # filename = 'PDB_resolution_method.pickle'
        filename = pdb + '_protein.pickle'

        self.logger.info("filename: filename: %s" % filename)
        return os.path.join("/usr/local/pipeline/hub-core/pickle-FR3D",filename)

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,
                            mod.UnitInfo.pdb_id,
                            mod.UnitInfo.number,
                            mod.UnitInfo.unit_type_id).\
                            filter(mod.UnitInfo.unit_type_id == 'aa').distinct()
        pdb_ids = []
        for row in query:
            pdb_ids.append(row.pdb_id)
        return pdb_ids[0:500]


    def data(self, pdb, **kwargs):
        """
            Look up all the existed pdbs to process.  Ignores the pdb input.
        """
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,
                            mod.UnitInfo.pdb_id,
                            mod.UnitInfo.number,
                            mod.UnitInfo.unit_type_id,
                            mod.UnitCenters.name,
                            mod.UnitCenters.x,
                            mod.UnitCenters.y,
                            mod.UnitCenters.z).\
                            join(mod.UnitCenters, mod.UnitInfo.unit_id == mod.UnitCenters.unit_id).\
                            filter(mod.UnitInfo.pdb_id==pdb).\
                            filter(or_(mod.UnitCenters.name == 'aa_fg',mod.UnitCenters.name == 'aa_backbone'))
        ids = []
        chainIndice = []
        aa_fg = []
        aa_backbone = []
        for row in query:
            # result['ids'].append(row.unit_id)
            # result['chainIndices'].append(str(row.number).replace('L',''))
            # result['centers'].append([row.x,row.y,row.z])
            if row.unit_id not in ids:
                ids.append(row.unit_id)
                chainIndice.append(int(str(row.number).replace('L','')))
                if row.name == 'aa_fg':
                    aa_fg.append((np.asarray([row.x,row.y,row.z])))
                    aa_backbone.append(np.asarray([np.NaN,np.NaN,np.NaN]))
                elif row.name == 'aa_backbone':
                    aa_backbone.append((np.asarray([row.x,row.y,row.z])))
                    aa_fg.append(np.asarray([np.NaN,np.NaN,np.NaN]))
            else:
                if row.name == 'aa_fg':
                    aa_fg[len(aa_fg)-1]=(np.asarray([row.x,row.y,row.z]))
                elif row.name == 'aa_backbone':
                    aa_backbone[len(aa_fg)-1]=(np.asarray([row.x,row.y,row.z]))





        # print(ids[:10],chainIndice[:10],aa_fg[:10],aa_backbone[:10])

        return [ids,chainIndice,aa_fg,aa_backbone]



    def process(self, pdb, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        webroot = self.config['locations']['fr3d_pickle_base'] + "/units/"

        filename = self.filename(pdb)

        pinfo = self.data(pdb)


        with open(filename, 'wb') as fh:
            self.logger.info("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

        os.system("rsync -u %s %s" % (filename, webroot))

        #copy to unit folder 
        copy_file = os.path.join("/var/www/html/units/",pdb + '_protein.pickle')
        with open(copy_file, 'wb') as fh:
            self.logger.info("process: copy_file open: %s" % copy_file)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

        os.system("rsync -u %s %s" % (copy_file, webroot))


        pass


