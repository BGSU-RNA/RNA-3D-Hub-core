"""
Module for export of pairs data in pickle format for FR3D.
"""

import numpy as np
import os
import pickle

from pymotifs import core
from pymotifs import models as mod

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.ife.info import Loader as IfeInfoLoader
from pymotifs.interactions.pairwise import Loader as InteractionLoader
# from pymotifs.interactions.annotate_python import Loader as AnnotationLoader

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
    # dependencies = set([ChainLoader, InteractionLoader, AnnotationLoader, IfeInfoLoader])
    dependencies = set([ChainLoader, InteractionLoader, IfeInfoLoader])

    def has_data(self, pdb, *args, **kwargs):
        self.logger.debug("has_data: pdb: %s" % str(pdb))
        filename = self.filename(pdb)

        #self.logger.info("Replacing %s" % filename)
        #return False

        if os.path.exists(filename):
            self.logger.info("Filename %s exists" % filename)
            return True
        else:
            self.logger.info("Filename %s is missing" % filename)
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

        # TO DO: change RNA to NA

        filename = pdb + '_RNA_pairs.pickle'

        #self.logger.info("filename: %s" % filename)

        return os.path.join("pickle-FR3D",filename)


    def data(self, pdb, **kwargs):
        """
        Get all pairs for the given PDB and format them for
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
            self.logger.info("data: building UnitPairs query")

            iui1 = aliased(mod.UnitInfo)
            iui2 = aliased(mod.UnitInfo)
            iupf = mod.UnitPairsFlanking
            iupi = mod.UnitPairsInteractions

            subqueryI = session.query(iupi.unit_id_1,
                               iupi.unit_id_2,
                               iui1.model,
                               iui1.chain,
                               iui1.chain_index,
                               iui2.model,
                               iui2.chain,
                               iui2.chain_index,
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
                    filter(iupi.program == 'matlab').\
                    filter(iui1.unit_type_id == 'rna').\
                    filter(iui2.unit_type_id == 'rna').\
                    filter(iupi.pdb_id == pdb)

            fui1 = aliased(mod.UnitInfo)
            fui2 = aliased(mod.UnitInfo)
            fupf = mod.UnitPairsFlanking
            fupi = mod.UnitPairsInteractions

            subqueryF = session.query(fupf.unit_id_1,
                               fupf.unit_id_2,
                               fui1.model,
                               fui1.chain,
                               fui1.chain_index,
                               fui2.model,
                               fui2.chain,
                               fui2.chain_index,
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
                    filter(fupi.program == 'matlab').\
                    filter(fui1.unit_type_id == 'rna').\
                    filter(fui2.unit_type_id == 'rna').\
                    filter(fupf.pdb_id == pdb)

            query = subqueryI.union(subqueryF).order_by(iui1.model, iui1.chain, iui1.chain_index,
                               iui2.model, iui2.chain, iui2.chain_index)

            self.logger.debug("data: query built: %s" % str(query))

            count = query.count()
            if not count:
                self.logger.warning("No interactions found for %s", pdb)
            else:
                self.logger.info("Found %s interactions for %s", count, pdb)

            # TO DO:  implement parsing logic from pairFileParsing.py here
            # and return that data instead of the full load.

            interactionToPair = defaultdict(list)

            #constraintList = ['f_lwbp', 'f_stacks', 'f_bphs', 'f_brbs']

            for result in query:
                uid1 = result.unit_id_1
                uid2 = result.unit_id_2

                if uid1 == 'placeholder':
                    continue

                # FUTURE: iterate over a list of column names?

                if( result.f_lwbp is not None and len(result.f_lwbp) > 2):
                    interactionToPair[result.f_lwbp].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_lwbp, uid1, uid2, result.f_crossing))

                if( result.f_stacks is not None and len(result.f_stacks) > 2):
                    interactionToPair[result.f_stacks].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_stacks, uid1, uid2, result.f_crossing))

                if( result.f_bphs is not None and len(result.f_bphs) > 2):
                    interactionToPair[result.f_bphs].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_bphs, uid1, uid2, result.f_crossing))

                if( result.f_brbs is not None and len(result.f_brbs) > 2):
                    interactionToPair[result.f_brbs].append((uid1, uid2, result.f_crossing))
                    self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_brbs, uid1, uid2, result.f_crossing))

                #for constraint in constraintList:
                #    self.logger.info("constraint: %s" % constraint)
                #    #self.logger.info("constraint value: %s" % result.constraint.value())
                #    #if len(result.constraint.value()) > 2:
                #    #    interactionToPair[constraint].append((uid1, uid2, result.f_crossing))
                #    #    self.logger.info("type: units/constraint: %s : %s, %s / %s" % (result.getAttr(constraint), uid1, uid2, result.f_crossing))
                #    pass

                if result.flanking == "1":
                    interactionToPair["bSS"].append((uid1, uid2, None))
                    self.logger.debug("type: units/constraint: bSS : %s, %s / None" % (uid1, uid2))

                pass


            # categories to get from PairAnnotations table
            categories = ['sO']

            iann = mod.PairAnnotations

            query3 = session.query(iann.unit_id_1,
                               iann.unit_id_2,
                               iann.annotation,
                               iann.crossing).\
                    filter(iann.category.in_(categories)).\
                    filter(iann.pdb_id == pdb)

            for result in query3:
                uid1 = result.unit_id_1
                uid2 = result.unit_id_2

                if (result.annotation is not None and len(result.annotation) > 2):
                    interactionToPair[result.annotation].append((uid1, uid2, result.crossing))
                    #self.logger.info("annotation: %s %s %s %s" % (uid1, result.annotation, uid2, result.crossing))

            return interactionToPair


    def process(self, entry, **kwargs):
        """
        Load pairwise interaction data for the given PDB file.

        Parameters
        ----------
        entry : object
            The entry to process, should be a PDB file.
        **kwargs : dict
            Generic keyword arguments.
        """

        pinfo = self.data(entry)

        self.logger.debug("process: raw data: %s" % pinfo)

        filename = self.filename(entry)

        with open(filename, 'wb') as fh:
            self.logger.info("writing pickle file: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

        # public directory where WebFR3D will access the .pickle files
        webroot = self.config['locations']['fr3d_pickle_base'] + "/pairs/"

        os.system("rsync -u %s %s" % (filename, webroot))
        self.logger.info("rsync -u %s %s" % (filename, webroot))