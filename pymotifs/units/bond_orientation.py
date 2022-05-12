"""A module to compute and store glycosidic bond orientations in the database.
"""

import numpy as np
import math

from pymotifs import core
from pymotifs import models as mod
from pymotifs.units.info import Loader as InfoLoader

from fr3d.modified_parent_mapping import modified_nucleotides
from NA_unit_annotation import annotate_bond_orientation

class Loader(core.SimpleLoader):
    """A class to compute glycosidic bond orientations and load them into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True

    def query(self, session, pdb):
        """Create a query to look up glycosidic bond orientations

        :session: The session object to use.
        :pdb: The pdb id to query for.
        :returns: A query to get bond orientations.
        """

        return session.query(mod.UnitOrientations).\
            filter_by(pdb_id=pdb)


    def to_process(self, pdbs, **kwargs):
        """
        Find pdb files with no bond orientations computed.

        Parameters
        ----------
        pdbs : list of strings

        Returns
        -------
        pdbs_to_compute : list of strings, can be empty

        """

        # when just one file to process, pdbs is a string
        if isinstance(pdbs,str):
            pdbs = [pdbs]

        with self.session() as session:
            query = session.query(mod.UnitOrientations.pdb_id).\
                   distinct()

            pdbs_computed = set([r.pdb_id for r in query])

        pdbs_to_compute = sorted(set(pdbs) - pdbs_computed)

        # make sure to return at least one file name, o/w dispatcher complains
        if len(pdbs_to_compute) == 0:
            pdbs_to_compute = [pdbs[0]]

        self.logger.info("Found %d files to process" % len(pdbs_to_compute))

        return pdbs_to_compute

    def data(self, pdb, **kwargs):
        """
        Load the structure, compute the orientations, save.

        :pdb: The pdb to process
        :yields: Yields angle chi_degree and an annotation
        """

        # load the 3D structure file
        structure = self.structure(pdb)

        bond_orientations = annotate_bond_orientation(structure,pdb)

        for b in bond_orientations:
            yield mod.UnitOrientations(unit_id=b["unit_id"],
                                    pdb_id=b["pdb_id"],
                                    orientation=b["orientation"],
                                    chi_degree=b["chi_degree"])
