"""Module for export of pairs data in pickle format for FR3D.
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

from collections import defaultdict

from os import path

from sqlalchemy import and_
from sqlalchemy.orm import aliased
from sqlalchemy.sql import select
from sqlalchemy.sql import union


class Exporter(core.Loader):
    """Export pairs data in pickle format, one file per PDB.
    """


    # General Setup
    compressed = False 
    mark = False 
    dependencies = set([ChainLoader, CentersLoader, RotationsLoader, 
                        PositionLoader, IfeInfoLoader, MappingLoader])


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


    def filename(self, pdb, **kwargs):
        """Create the filename for the given PDB.

        Parameters
        ----------
        pdb : string
            The PDB ID for which to create a file.

        Returns
        -------
        filename : str
            The path to write to.
        """

        # TO DO: put the important directories into the config

        filename = pdb + '_RNA_pairs.pickle'

        self.logger.info("filename: filename: %s" % filename)

        return os.path.join("pickle-FR3D",filename)


    def data(self, pdb, **kwargs):
        """Get all pairs for the given PDB and format them for
        convenient use by FR3D.

        Parameters
        ----------
        pdb : string
            The PDB ID for which to look up pairs data.

        Returns
        -------
        resultset : list of dicts
        """

        with self.session() as session:
            self.logger.info("data: Inside data retrieval routine")
            self.logger.info("data: building query")

            fui1 = aliased(mod.UnitInfo)
            fui2 = aliased(mod.UnitInfo)
            fupf = mod.UnitPairsFlanking
            fupi = mod.UnitPairsInteractions

            iui1 = aliased(mod.UnitInfo)
            iui2 = aliased(mod.UnitInfo)
            iupf = mod.UnitPairsFlanking
            iupi = mod.UnitPairsInteractions

            subqueryI = session.query(iupi.unit_id_1,
                               iupi.unit_id_2,
                               iupi.f_lwbp,
                               iupi.f_stacks,
                               iupi.f_bphs,
                               iupi.f_brbs,
                               iupi.f_crossing,
                               iupf.flanking).\
                    join(iui1, iupi.unit_id_1 == iui1.unit_id).\
                    join(iui2, iupi.unit_id_2 == iui2.unit_id).\
                    outerjoin(iupf, and_(iupi.unit_id_1 == iupf.unit_id_1, iupi.unit_id_2 == iupf.unit_id_2, 
                              iupi.pdb_id == iupf.pdb_id)).\
                    filter(iupi.unit_id_1.isnot(None)).\
                    filter(iupi.unit_id_2.isnot(None)).\
                    filter(iui1.unit_type_id == 'rna').\
                    filter(iui2.unit_type_id == 'rna').\
                    filter(iupi.pdb_id == pdb)

            subqueryF = session.query(fupf.unit_id_1,
                               fupf.unit_id_2,
                               fupi.f_lwbp,
                               fupi.f_stacks,
                               fupi.f_bphs,
                               fupi.f_brbs,
                               fupi.f_crossing,
                               fupf.flanking).\
                    join(fui1, fupf.unit_id_1 == fui1.unit_id).\
                    join(fui2, fupf.unit_id_2 == fui2.unit_id).\
                    outerjoin(fupi, and_(fupi.unit_id_1 == fupf.unit_id_1, fupi.unit_id_2 == fupf.unit_id_2, 
                              fupi.pdb_id == fupf.pdb_id)).\
                    filter(fupf.unit_id_1.isnot(None)).\
                    filter(fupf.unit_id_2.isnot(None)).\
                    filter(fui1.unit_type_id == 'rna').\
                    filter(fui2.unit_type_id == 'rna').\
                    filter(fupf.pdb_id == pdb)

            query = subqueryI.union(subqueryF)

            self.logger.debug("data: query built: %s" % str(query))

            count = query.count()
            if not count:
                self.logger.warning("No interactions found for %s", pdb)
            else:
                self.logger.info("Found %s interactions for %s", count, pdb)

            # TO DO:  implement parsing logic from pairFileParsing.py here
            # and return that data instead of the full load.

            constraintToPair = defaultdict(list)

            #constraintList = ['f_lwbp', 'f_stacks', 'f_bphs', 'f_brbs']

            for result in query:
                uid1 = result.unit_id_1
                uid2 = result.unit_id_2

                # FUTURE: iterate over a list of column names?

                if( result.f_lwbp is not None and len(result.f_lwbp) > 2):
                    constraintToPair[result.f_lwbp].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_lwbp, uid1, uid2, result.f_crossing))

                if( result.f_stacks is not None and len(result.f_stacks) > 2):
                    constraintToPair[result.f_stacks].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_stacks, uid1, uid2, result.f_crossing))

                if( result.f_bphs is not None and len(result.f_bphs) > 2):
                    constraintToPair[result.f_bphs].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_bphs, uid1, uid2, result.f_crossing))

                if( result.f_brbs is not None and len(result.f_brbs) > 2):
                    constraintToPair[result.f_brbs].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_brbs, uid1, uid2, result.f_crossing))

                #for constraint in constraintList:
                #    self.logger.info("constraint: %s" % constraint)
                #    #self.logger.info("constraint value: %s" % result.constraint.value())
                #    #if len(result.constraint.value()) > 2:
                #    #    constraintToPair[constraint].append((uid1, uid2, result.f_crossing))
                #    #    self.logger.info("type: units/constraint: %s : %s, %s / %s" % (result.getAttr(constraint), uid1, uid2, result.f_crossing))
                #    pass

                if result.flanking == "1":
                    constraintToPair["bSS"].append((uid1, uid2, None))
                    self.logger.debug("type: units/constraint: bSS : %s, %s / None" % (uid1, uid2))

                pass

            #return [row2dict(result) for result in query] # full load data for testing
            return constraintToPair


    #def to_process(self, pdbs, **kwargs):
    #    """Process the input list of PDBs.

    #    Parameters
    #    ----------
    #    pdbs : list
    #        To process.

    #    Returns
    #    -------
    #    (pdb_id, model, chain) : tuple
    #        The components of the IFE-chains to be processed.
    #    """

    #    #with self.session() as session:
    #    #    query = session.query(
    #    #               mod.UnitInfo.pdb_id,
    #    #               mod.UnitInfo.model,
    #    #               mod.UnitInfo.chain
    #    #           ).\
    #    #           distinct().\
    #    #           filter(mod.UnitInfo.unit_type_id == 'rna')

    #    #    return [(r.pdb_id, r.model, r.chain) for r in query]

    #    #pass
    #    return pdbs


    def process(self, entry, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        entry : object
            The entry to process.
        **kwargs : dict
            Generic keyword arguments.
        """

        filename = self.filename(entry)

        pinfo = self.data(entry)

        self.logger.debug("process: raw data: %s" % pinfo)

        with open(filename, 'wb') as fh:
            self.logger.info("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

        pass

