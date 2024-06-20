"""
Module for export of unit annotations like syn/anti
in pickle format for FR3D.
"""

import numpy as np
import os
import pickle

from pymotifs import core
from pymotifs import models as mod

from pymotifs.units.bond_orientation import Loader as OrientationLoader

from os import path


class Exporter(core.Loader):
    """Export unit data in pickle format, one file per
    IFE-chain.
    """

    # General Setup
    compressed = False
    mark = False
    dependencies = set([OrientationLoader])

    def has_data(self, entry, *args, **kwargs):
        #self.logger.info("has_data: entry: %s" % str(entry))
        filename = self.filename(entry)
        #self.logger.info("has_data: filename: %s" % filename)

        # force re-writing of .pickle files
        #return False

        if os.path.exists(filename) is True:
            self.logger.info("has_data: filename %s exists" % filename)
            return True
        else:
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

        return os.path.join("pickle-FR3D",chain_string + "_NA_unit_annotations.pickle")


    def to_process(self, pdbs, **kwargs):
        """
        Look up the list of chains to process.

        Parameters
        ----------
        pdbs : list

        Returns
        -------
        (pdb_id, model, chain) : tuple
            The components of the IFE-chains to be processed.
        """

        molecule_types = ['rna','dna']

        # speed up the query by limiting to pdb_id in pdbs
        with self.session() as session:
            query = session.query(
                       mod.UnitInfo.pdb_id,
                       mod.UnitInfo.model,
                       mod.UnitInfo.chain).\
                   distinct().\
                   filter(mod.UnitInfo.unit_type_id.in_(molecule_types)).\
                   filter(mod.UnitInfo.pdb_id.in_(pdbs))

            return [(r.pdb_id, str(r.model), r.chain) for r in query if r.pdb_id in pdbs]


    def data(self, ichain, **kwargs):
        """
        Get all unit annotations for the given IFE-chain
        and format them for convenient use by FR3D.

        Parameters
        ----------
        ichain : tuple
            The IFE-chain for which to look up centers and rotation data.

        Returns
        -------
        resultset : list of dicts
        """

        pdb = ichain[0]
        mdl = str(ichain[1])
        chn = ichain[2]

        unit_id_to_annotations = {}

        with self.session() as session:
            query = session.query(mod.UnitAnnotations.unit_id,
                               mod.UnitAnnotations.category,
                               mod.UnitAnnotations.value).\
                     filter(mod.UnitAnnotations.pdb_id == pdb)

            # loop over lines from the table
            for row in query:
                fields = row.unit_id.split("|")
                # match model and chain
                if fields[1] == mdl and fields[2] == chn:
                    if row.unit_id in unit_id_to_annotations:
                        unit_id_to_annotations[row.unit_id][row.category] = row.value
                    else:
                        unit_id_to_annotations[row.unit_id] = {}
                        unit_id_to_annotations[row.unit_id][row.category] = row.value

        return unit_id_to_annotations


    def process(self, entry, **kwargs):
        """Load unit annotation data for the given IFE-chain.

        Parameters
        ----------
        entry : object
            The entry to process.
        **kwargs : dict
            Generic keyword arguments.
        """

        # "local" location of file to write out
        filename = self.filename(entry)

        annotations = self.data(entry)

        with open(filename, 'wb') as fh:
            self.logger.debug("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(annotations, fh, 2)

        # location of files made available via the web server
        webroot = self.config['locations']['fr3d_pickle_base'] + "/units/"

        # synchronize current file between hub-core and web server directories
        os.system("rsync -u %s %s" % (filename, webroot))
        self.logger.debug("rsync -u %s %s" % (filename, webroot))

