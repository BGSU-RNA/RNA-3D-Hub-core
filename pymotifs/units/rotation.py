"""A module to load unit rotation matrices for RNA into the database.
"""

import pymotifs.core as core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitRotations


class Loader(core.SimpleLoader):
    """A class to load rotation matrices into the database. This will load the
    RNA rotation matrices into the database.
    """

    def query(self, session, pdb):
        """Create a query to lookup the rotation matrices.

        :session: The session object to use.
        :pdb: The pdb id to query for.
        :returns: A query to get rotation matrices.
        """

        return session.query(UnitRotations).\
            join(UnitInfo, UnitInfo.id == UnitRotations.id).\
            filter(UnitInfo.pdb == pdb)

    def data(self, pdb):
        """Get the rotation matrices for all RNA residues in the given pdb.

        :pdb: The pdb to process
        :yields: Yields a series of rotation matrices.
        """

        structure = self.structure(pdb)
        for residue in structure.residues():
            if hasattr(residue, 'rotation_matrix'):
                matrix = residue.rotation_matrix
                yield UnitRotations(id=residue.unit_id(),
                                    r1c1=matrix[0, 0],
                                    r1c2=matrix[0, 1],
                                    r1c3=matrix[0, 2],
                                    r2c1=matrix[1, 0],
                                    r2c2=matrix[1, 1],
                                    r2c3=matrix[1, 2],
                                    r3c1=matrix[2, 0],
                                    r3c2=matrix[2, 1],
                                    r3c3=matrix[2, 2])
