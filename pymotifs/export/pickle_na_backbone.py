"""Module for export of nt backbone centers in pickle format for FR3D.
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
from sqlalchemy.orm import aliased

from os import path


class Exporter(core.Loader):
    """Export unit data in pickle format, one file per
    IFE-chain.
    """

    # General Setup
    compressed = False
    mark = False
    dependencies = set([ChainLoader, CentersLoader,
                        PositionLoader, IfeInfoLoader, MappingLoader])


    def has_data(self, entry, *args, **kwargs):
        self.logger.info("has_data: entry: %s" % str(entry))
        filename = self.filename(entry)
        self.logger.info("has_data: filename: %s" % filename)
        if os.path.exists(filename) is True:
            self.logger.info("has_data: filename %s exists" % filename)
            return True
        self.logger.info("has_data: filename %s is missing" % filename)
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

        # write directly to web directory, don't keep two copies of the file
        return os.path.join(self.config['locations']['fr3d_pickle_base'],"units",chain_string + "_NA_phosphate_sugar.pickle")


    def data(self, ichain, **kwargs):
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
            self.logger.debug("na_backbone: Inside data retrieval routine")
            self.logger.debug("na_backbone: building query")

            unit_type_ids = ['rna','dna']
            centers = ['nt_phosphate','nt_sugar']

            centers1 = aliased(mod.UnitCenters)   # first reference for nt_phosphate
            centers2 = aliased(mod.UnitCenters)   # second reference for nt_sugar

            query = session.query(mod.UnitInfo.unit_id,
                               mod.ExpSeqPosition.index.label('position_order'),
                               centers1.x,
                               centers1.y,
                               centers1.z,
                               centers2.x.label('xx'),
                               centers2.y.label('yy'),
                               centers2.z.label('zz')).\
                     distinct().\
                     join(centers1, mod.UnitInfo.unit_id == centers1.unit_id).\
                     join(centers2, mod.UnitInfo.unit_id == centers2.unit_id).\
                     join(mod.ExpSeqUnitMapping, mod.UnitInfo.unit_id == mod.ExpSeqUnitMapping.unit_id).\
                     join(mod.ExpSeqPosition, mod.ExpSeqUnitMapping.exp_seq_position_id == mod.ExpSeqPosition.exp_seq_position_id).\
                     filter(mod.UnitInfo.unit_type_id == "rna").\
                     filter(centers1.name == 'nt_phosphate').\
                     filter(centers2.name == 'nt_sugar').\
                     filter(mod.UnitInfo.pdb_id == pdb).\
                     filter(mod.UnitInfo.model == str(mdl)).\
                     filter(mod.UnitInfo.chain == chn).\
                     order_by(mod.ExpSeqPosition.index)

            self.logger.debug("na_backbone: query built")

            units = []
            order = []
            phos  = []
            sugar = []
            rttns = []
            rsset = []

            for row in query:
                units.append(row.unit_id)
                order.append(row.position_order)
                phos.append( np.asarray([row.x, row.y, row.z]))
                sugar.append(np.asarray([row.xx, row.yy, row.zz]))

            rsset = [ units, order, phos, sugar ]

            self.logger.debug("na_backbone: units: %s" % str(units))
            self.logger.debug("na_backbone: order: %s" % str(order))
            self.logger.debug("na_backbone: phos : %s" % str(phos))
            self.logger.debug("na_backbone: sugar: %s" % str(sugar))
            self.logger.debug("na_backbone: rsset: %s" % str(rsset))

            return rsset


    def to_process(self, pdbs, **kwargs):
        """Look up the list of IFE-chains to process.  Ignores the pdbs input.

        Parameters
        ----------
        pdbs : list
            Ignored unless it's a short list

        Returns
        -------
        (pdb_id, model, chain) : tuple
            The components of the IFE-chains to be processed.
        """

        self.logger.info("Given %d pdbs to process" % len(pdbs))

        with self.session() as session:
            query = session.query(
                       mod.UnitInfo.pdb_id,
                       mod.UnitInfo.model,
                       mod.UnitInfo.chain
                   ).\
                   distinct().\
                   filter(mod.UnitInfo.unit_type_id == 'rna')

            # if only a few pdbs are requested, focus on those, otherwise do all
            if len(pdbs) < 10:
                all_pdbs = [(r.pdb_id, r.model, r.chain) for r in query if r.pdb_id in pdbs]
            else:
                all_pdbs = [(r.pdb_id, r.model, r.chain) for r in query]

            self.logger.info("Found %d pdbs to process" % len(all_pdbs))

            return all_pdbs


    def process(self, entry, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        entry : object
            The entry to process.
        **kwargs : dict
            Generic keyword arguments.
        """

        webroot = self.config['locations']['fr3d_pickle_base'] + "/units/"

        filename = self.filename(entry)

        uinfo = self.data(entry)

        with open(filename, 'wb') as fh:
            self.logger.debug("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(uinfo, fh, 2)

#        os.system("rsync -u %s %s" % (filename, webroot))
#        self.logger.debug("rsync -u %s %s" % (filename, webroot))

