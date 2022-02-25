"""Module for export of resolution and methods of existing pdb files
"""

import numpy as np
import os
import pickle

from pymotifs import core
from pymotifs import models as mod

# from pymotifs.chains.info import Loader as ChainLoader
# from pymotifs.units.centers import Loader as CentersLoader
# from pymotifs.units.rotation import Loader as RotationsLoader
# from pymotifs.exp_seq.mapping import Loader as MappingLoader
# from pymotifs.exp_seq.positions import Loader as PositionLoader
# from pymotifs.ife.info import Loader as IfeInfoLoader
from pymotifs.pdbs.info import Loader as InfoLoader

from collections import defaultdict

from os import path

from sqlalchemy import and_
from sqlalchemy.orm import aliased
from sqlalchemy.sql import select
from sqlalchemy.sql import union


class Exporter(core.Loader):
    """Export pairs data in pickle format.
    """


    # General Setup
    compressed = False 
    mark = False 
    # dependencies = set([InfoLoader])
    dependencies = set()


    def has_data(self, *args, **kwargs):
        filename = self.filename()
        if os.path.exists(filename) is True:
            return True
        return False


    def remove():
        pass


    def filename(self):
        """Create the filename for the given PDB.

        Parameters
        ----------
        pdb : 

        Returns
        -------
        filename : str
            The path to write to.
        """

        # TO DO: put the important directories into the config

        filename = 'PDB_resolution_method.pickle'

        self.logger.info("filename: filename: %s" % filename)

        return os.path.join("logs",filename)


    def data(self):
        """
            This function is a query function.
        """

        with self.session() as session:
            query = session(mod.PdbId.pdb_id,
                            mod.PdbId.resolution,
                            mod.PdbId.experimental_technique)
        # result = {row.pdb_id:{'resolution': row.resolution,'method': row.experimental_technique } for row in query}

        result = {}
        for row in query:
            result[row.pdb_id] = {'resolution': row.resolution,'method': row.experimental_technique }
        return result



    def process(self, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        # webroot = self.config['locations']['fr3d_pickle_base'] + "/pairs/"

        filename = self.filename()

        pinfo = self.data()


        with open(filename, 'wb') as fh:
            self.logger.info("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

        # os.system("rsync -u %s %s" % (filename, webroot))


        pass

