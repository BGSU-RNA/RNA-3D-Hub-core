"""A module to load unit rotation matrices for RNA into the database.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    """A class to load rotation matrices into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True

    def query(self, session, pdb):
        """Create a query to lookup the rotation matrices.

        :session: The session object to use.
        :pdb: The pdb id to query for.
        :returns: A query to get rotation matrices.
        """

        return session.query(mod.UnitRotations).\
            filter_by(pdb_id=pdb)

    def data(self, pdb, **kwargs):
        """Get the rotation matrices for all RNA residues in the given pdb.

        :pdb: The pdb to process
        :yields: Yields a series of rotation matrices.
        """

        structure = self.structure(pdb)

        # running infer_hydrogens, but now that is done automatically
#        structure.infer_hydrogens()

        for residue in structure.residues():

            if hasattr(residue, 'rotation_matrix'):
                matrix = residue.rotation_matrix
                # if there are not enough atoms for the rotation matrix, it will be None
                if matrix is not None:
                    yield mod.UnitRotations(unit_id=residue.unit_id(),
                                            pdb_id=pdb,
                                            cell_0_0=float(matrix[0, 0]),
                                            cell_0_1=float(matrix[0, 1]),
                                            cell_0_2=float(matrix[0, 2]),
                                            cell_1_0=float(matrix[1, 0]),
                                            cell_1_1=float(matrix[1, 1]),
                                            cell_1_2=float(matrix[1, 2]),
                                            cell_2_0=float(matrix[2, 0]),
                                            cell_2_1=float(matrix[2, 1]),
                                            cell_2_2=float(matrix[2, 2]))
