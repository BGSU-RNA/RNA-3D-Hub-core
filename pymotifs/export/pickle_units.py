"""Module for export of unit center/rotation data 
in pickle format for FR3D.
"""

import numpy as np
import os
import pickle

from pymotifs import core
from pymotifs import models as mod

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.units.centers import Loader as CentersLoader
from pymotifs.units.rotation import Loader as RotationsLoader
from pymotifs.exp_seq.mapping import Loader as MappingLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader
from pymotifs.ife.info import Loader as IfeInfoLoader


class Exporter(core.Loader):
    """Export unit data in pickle format, one file per 
    IFE-chain.
    """


    # General Setup
    compressed = False 
    mark = False 
    dependencies = set([ChainLoader, CentersLoader, RotationsLoader, 
                        PositionLoader, IfeInfoLoader, MappingLoader])


    def has_data(self, *args, **kwargs):
        return False


    def remove():
        pass


    def filename(self, ichain, **kwargs):
        """Create the filename for the given IFE-chain.

        Parameters
        ----------
        ichain : tuple
            The components (PDB ID, model, and chain) of the 
            IFE-chain name (hyphen-deliminted concatenation) 
            for which to create a file.

        Returns
        -------
        filename : str
            The path to write to.
        """

        # TO DO: put the important directories into the config

        pdb = ichain[0]
        mdl = ichain[1]
        chn = ichain[2]

        chain_string = pdb + '-' + str(mdl) + '-' + chn

        self.logger.debug("filename: chain_string: %s" % chain_string)

        return os.path.join("pickle-FR3D",chain_string + "_RNA.pickle")


    def cenrot(self, ichain):
        """Get all unit listings for the given IFE-chain, centers and
        rotations, and format them for convenient use by FR3D.

        Parameters
        ----------
        ichain : tuple
            The IFE-chain for which to look up centers and rotation data.

        Returns
        -------
        resultset : list of dicts
        """

        pdb = ichain[0]
        mdl = ichain[1]
        chn = ichain[2]

        with self.session() as session:
            self.logger.debug("cenrot: Inside data retrieval routine")
            self.logger.debug("cenrot: building query")

            query = session.query(mod.UnitInfo.unit_id,
                               mod.ExpSeqPosition.index.label('position_order'),
                               mod.UnitCenters.x,
                               mod.UnitCenters.y,
                               mod.UnitCenters.z,
                               mod.UnitRotations.cell_0_0,
                               mod.UnitRotations.cell_0_1,
                               mod.UnitRotations.cell_0_2,
                               mod.UnitRotations.cell_1_0,
                               mod.UnitRotations.cell_1_1,
                               mod.UnitRotations.cell_1_2,
                               mod.UnitRotations.cell_2_0,
                               mod.UnitRotations.cell_2_1,
                               mod.UnitRotations.cell_2_2).\
                     distinct().\
                     join(mod.UnitCenters, mod.UnitInfo.unit_id == mod.UnitCenters.unit_id).\
                     join(mod.UnitRotations, mod.UnitInfo.unit_id == mod.UnitRotations.unit_id).\
                     join(mod.ExpSeqUnitMapping, mod.UnitInfo.unit_id == mod.ExpSeqUnitMapping.unit_id).\
                     join(mod.ExpSeqPosition, mod.ExpSeqUnitMapping.exp_seq_position_id == mod.ExpSeqPosition.exp_seq_position_id).\
                     filter(mod.UnitInfo.unit_type_id == 'rna').\
                     filter(mod.UnitCenters.name == 'base').\
                     filter(mod.UnitInfo.pdb_id == pdb).\
                     filter(mod.UnitInfo.model == str(mdl)).\
                     filter(mod.UnitInfo.chain == chn).\
                     order_by(mod.ExpSeqPosition.index)

            self.logger.debug("cenrot: query built")

            units = []
            order = []
            cntrs = []
            rttns = [] 
            rsset = []

            for row in query:
                units.append(row.unit_id)
                order.append(row.position_order)
                cntrs.append(np.asarray([row.x, row.y, row.z]))
                rttns.append(np.asarray([np.asarray([row.cell_0_0, row.cell_0_1, row.cell_0_2]),
                             np.asarray([row.cell_1_0, row.cell_1_1, row.cell_1_2]),
                             np.asarray([row.cell_2_0, row.cell_2_1, row.cell_2_2])]))
                
            rsset = [ units, order, cntrs, rttns ]

            self.logger.debug("cenrot: units: %s" % str(units))
            self.logger.debug("cenrot: order: %s" % str(order))
            self.logger.debug("cenrot: cntrs: %s" % str(cntrs))
            self.logger.debug("cenrot: rttns: %s" % str(rttns))
            self.logger.debug("cenrot: rsset: %s" % str(rsset))

            return rsset


    def data(self, session, cenrot, filename, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        cenrot : list of dicts
            The centers/rotations data to write

        filename : text
            The name of the output pickle file.


        Returns
        -------
        pickle : pickle
            A pickle file containing data for all units in the input chain.
        """

        pass


    def to_process(self, pdbs, **kwargs):
        """Look up the list of IFE-chains to process.  Ignores the pdbs input.

        Parameters
        ----------
        pdbs : list
            Ignored.

        Returns
        -------
        (pdb_id, model, chain_name) : tuple
            The components of the IFE-chains to be processed.
        """

        with self.session() as session:
            query = session.query(
                       mod.IfeInfo.pdb_id,
                       mod.IfeInfo.model,
                       mod.ChainInfo.chain_name
                   ).\
                   distinct().\
                   join(mod.IfeChains,
                        mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
                   join(mod.ChainInfo,
                        mod.ChainInfo.chain_id == mod.IfeChains.chain_id).\
                   filter(mod.IfeInfo.model.isnot(None))

            return [(r.pdb_id, r.model, r.chain_name) for r in query]

        pass


    def process(self, entry, **kwargs):
        """Process this entry. In the case of loaders this will parse the data
        and put it into the database, exporters may go to the database and then
        generate the file. Inheriting classes must implement this.

        Parameters
        ----------
        entry : object
            The entry to process.
        **kwargs : dict
            Generic keyword arguments.
        """

        filename = self.filename(entry)

        uinfo = self.cenrot(entry)

        with open(filename, 'wb') as fh:
            self.logger.debug("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(uinfo, fh, 2)

        pass

