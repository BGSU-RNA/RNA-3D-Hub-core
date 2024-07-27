"""
Module for export of pairs data in pickle format for FR3D.
"""

from collections import defaultdict
import os
import pickle
from sqlalchemy import and_
from sqlalchemy.orm import aliased
from sqlalchemy.sql import select
from sqlalchemy.sql import union

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import DATA_FILE_DIRECTORY

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.ife.info import Loader as IfeInfoLoader
from pymotifs.interactions.loader import Loader as InteractionLoader


class Exporter(core.Loader):
    """
    Export pair data in pickle format, one file per PDB.
    """

    # General Setup
    compressed = False
    mark = False
    # dependencies = set([ChainLoader, InteractionLoader, AnnotationLoader, IfeInfoLoader])
    dependencies = set([ChainLoader, InteractionLoader, IfeInfoLoader])

    def has_data(self, pdb, *args, **kwargs):

        filename = self.filename(pdb)

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

        filename = pdb + '_NA_pairs.pickle'

        return os.path.join(DATA_FILE_DIRECTORY,'pairs',filename)


    def data(self, pdb, **kwargs):
        """
        Get all pairs for the given PDB and format them for
        convenient use by FR3D.
        Will include pairs involving RNA, DNA, hybrid, whatever is annotated.

        Parameters
        ----------
        pdb : string
            The PDB ID for which to look up pairs data.

        Returns
        -------
        resultset : list of dicts
        """

        with self.session() as session:
            iupi = mod.UnitPairsInteractions

            query = session.query(iupi.unit_id_1,
                               iupi.unit_id_2,
                               iupi.f_lwbp,
                               iupi.f_stacks,
                               iupi.f_bphs,
                               iupi.f_brbs,
                               iupi.f_crossing,
                               iupi.f_so,
                               iupi.f_coplanar,
                               iupi.f_sugar_ribose,
                               iupi.f_bss,
                               iupi.f_covalent).\
                    filter(iupi.unit_id_1.isnot(None)).\
                    filter(iupi.unit_id_2.isnot(None)).\
                    filter(iupi.program == 'python').\
                    filter(iupi.pdb_id == pdb)

            # count = query.count()
            # if not count:
            #     self.logger.info("No interactions found for %s", pdb)
            # else:
            #     self.logger.info("Found %s interactions for %s", count, pdb)

            interactionToPair = defaultdict(list)

            category_list = ['f_lwbp', 'f_stacks', 'f_bphs', 'f_brbs', 'f_so', 'f_coplanar', 'f_sugar_ribose', 'f_bss', 'f_covalent']
            for result in query:
                unit_id_1 = result.unit_id_1
                unit_id_2 = result.unit_id_2

                if unit_id_1 == 'placeholder':
                    continue

                for category in category_list:
                    interaction = getattr(result, category)
                    if interaction is not None and len(interaction) > 1:
                       interactionToPair[interaction].append((unit_id_1, unit_id_2, result.f_crossing))
                       self.logger.info("Recording: %s %s %s %s" % (unit_id_1, interaction, unit_id_2, result.f_crossing))

                # if result.f_lwbp is not None and len(result.f_lwbp) > 2:
                #     interactionToPair[result.f_lwbp].append((unit_id_1, unit_id_2, result.f_crossing))
                #     self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_lwbp, unit_id_1, unit_id_2, result.f_crossing))

                # if result.f_stacks is not None and len(result.f_stacks) > 2:
                #     interactionToPair[result.f_stacks].append((unit_id_1, unit_id_2, result.f_crossing))
                #     self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_stacks, unit_id_1, unit_id_2, result.f_crossing))

                # if result.f_bphs is not None and len(result.f_bphs) > 2:
                #     interactionToPair[result.f_bphs].append((unit_id_1, unit_id_2, result.f_crossing))
                #     self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_bphs, unit_id_1, unit_id_2, result.f_crossing))

                # if result.f_brbs is not None and len(result.f_brbs) > 2:
                #     interactionToPair[result.f_brbs].append((unit_id_1, unit_id_2, result.f_crossing))
                #     self.logger.debug("type: units/constraint: %s : %s, %s / %s" % (result.f_brbs, unit_id_1, unit_id_2, result.f_crossing))


            return interactionToPair


    def process(self, pdb, **kwargs):
        """
        Load pairwise interaction data for the given PDB file.

        Parameters
        ----------
        pdb : object
            The entry to process, should be a PDB file.
        **kwargs : dict
            Generic keyword arguments.
        """

        pinfo = self.data(pdb)

        filename = self.filename(pdb)

        with open(filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)
