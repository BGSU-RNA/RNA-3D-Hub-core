"""
Module for export of pairs data in pickle format for FR3D.
Takes about 17 minutes to export 20,000 files.
Files are stored in hub-core/data/pairs
Runs on just the current files in the PDB query.
"""

from collections import defaultdict
import os
import pickle

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
    dependencies = set([ChainLoader, InteractionLoader, IfeInfoLoader])


    # def to_process(self, pdbs, **kwargs):
    #     """
    #     Comment this out for weekly use!

    #     Return *all* PDB ids which have pairwise annotations.
    #     :param list pdbs: The list of pdb ids of interest on this run, ignored.
    #     :returns: A list of pdb ids to process.
    #     """

    #     with self.session() as session:
    #         query = session.query(mod.UnitPairsInteractions2024.pdb_id).\
    #         distinct()
    #         pdb_id_all = set([result.pdb_id for result in query])

    #     return sorted(pdb_id_all)


    def filename(self, pdb, **kwargs):
        """
        Create the filename for the given PDB.

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
            iupi = mod.UnitPairsInteractions2024

            query = session.query(iupi.unit_id_1,
                               iupi.unit_id_2,
                               iupi.f_lwbp_detail,
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
                    filter(iupi.pdb_id == pdb)

            interactionToPair = defaultdict(list)

            category_list = ['f_lwbp_detail', 'f_stacks', 'f_bphs', 'f_brbs', 'f_so', 'f_coplanar', 'f_sugar_ribose', 'f_bss', 'f_covalent']
            c = 0
            for result in query:
                unit_id_1 = result.unit_id_1
                unit_id_2 = result.unit_id_2

                if unit_id_1 == 'placeholder':
                    continue

                for category in category_list:
                    interaction = getattr(result, category)
                    if interaction is not None:
                        interactionToPair[interaction].append((unit_id_1, unit_id_2, result.f_crossing))
                        # self.logger.info("Recording: %s %s %s %s %s" % (category,unit_id_1, interaction, unit_id_2, result.f_crossing))
                        # print("Recording: %s %s %s %s %s" % (category,unit_id_1, interaction, unit_id_2, result.f_crossing))
                        c += 1

            self.logger.info("Recording %5d pairs from %s" % (c,pdb))

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
