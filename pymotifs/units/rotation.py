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
        structure.infer_hydrogens()
        for residue in structure.residues():
            if hasattr(residue, 'rotation_matrix'):
                matrix = residue.rotation_matrix
                yield mod.UnitRotations(unit_id=residue.unit_id(),
                                        pdb_id=pdb,
                                        cell_0_0=matrix[0, 0],
                                        cell_0_1=matrix[0, 1],
                                        cell_0_2=matrix[0, 2],
                                        cell_1_0=matrix[1, 0],
                                        cell_1_1=matrix[1, 1],
                                        cell_1_2=matrix[1, 2],
                                        cell_2_0=matrix[2, 0],
                                        cell_2_1=matrix[2, 1],
                                        cell_2_2=matrix[2, 2])
