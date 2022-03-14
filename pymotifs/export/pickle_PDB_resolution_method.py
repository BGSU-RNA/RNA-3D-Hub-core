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
    dependencies = set([InfoLoader])
    #dependencies = set()


    def has_data(self, *args, **kwargs):
        filename = self.filename()
        return False
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


        filename = 'PDB_resolution_method.pickle'

        self.logger.info("filename: filename: %s" % filename)

        return os.path.join("logs",filename)

    def to_process(self, pdbs, **kwargs):
        """Look up the list of IFE-chains to process.  Ignores the pdbs input.

        Parameters
        ----------
        pdbs : list
            Ignored.

        Returns
        -------
            all existed pdb ids.
        """
        return [0]
        with self.session() as session:
            query = session.query(
                       mod.PdbInfo.pdb_id
                   ).\
                   distinct()

            #return [r.pdb_id for r in query]



    def data(self, pdb, **kwargs):
        """
            This function is a query function.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.pdb_id,
                            mod.ChainInfo.chain_name,
                            mod.ChainInfo.entity_macromolecule_type,
                            mod.PdbInfo.resolution,
                            mod.PdbInfo.experimental_technique).\
                            join(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id)
        # result = {row.pdb_id:{'resolution': row.resolution,'method': row.experimental_technique } for row in query}
        print('am i here???????????????')
        result = {}
        result['chain'] = {}
        for row in query:
            #result[row.pdb_id] = {'resolution': row.resolution,'method': row.experimental_technique ,row.chain_name:row.entity_macromolecule_type}
            # result['pdb_id'] = row.pdb_id
            # result['resolution'] = row.resolution
            # result['method'] = row.experimental_technique
            # chaintype = row.entity_macromolecule_type
            # if chaintype in result['chain']:
            #     result['chain'][chaintype].append(row.chain_name)
            # else:
            #     result['chain'][chaintype] = [row.chain_name]
            # result['chain'].append()
            result[row.pdb_id][] = 
        print(result)

        return result



    def process(self, pdb, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        # webroot = self.config['locations']['fr3d_pickle_base'] + "/pairs/"

        filename = self.filename()

        pinfo = self.data(pdb)


        with open(filename, 'wb') as fh:
            self.logger.info("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

        # os.system("rsync -u %s %s" % (filename, webroot))


        pass

