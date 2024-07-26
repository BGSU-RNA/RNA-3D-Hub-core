"""
A module to load unit centers and rotation matrices for nucleotides into the database.
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    """
    A class to load unit centers and rotation matrices into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True

    def to_process(self, pdbs, **kwargs):
        """
        Do one query to determine which pdb ids do not have centers/rotations
        Since you need a center to have a rotation, just check for centers
        """

        with self.session() as session:
            query = session.query(mod.UnitCenters.pdb_id).distinct()

        existing_ids = set()
        for result in query:
            existing_ids.add(result.pdb_id)

        pdbs_to_check = set(pdbs) - existing_ids

        self.logger.info("Found %d existing ids and %d pdbs to check" % (len(existing_ids),len(pdbs_to_check)))

        if len(pdbs_to_check) == 0:
            raise core.Skip("All pdb ids already represented in unit_rotations table")

        return sorted(pdbs_to_check)

    def query(self, session, pdb):
        """
        Create a query to look up the unit centers.

        :session: The session object to use.
        :pdb: The pdb id to query for.
        :returns: A query to get rotation matrices.
        """

        return session.query(mod.UnitCenters).\
            filter_by(pdb_id=pdb)

    def data(self, pdb, **kwargs):
        """
        Load the structure, extract centers and rotation matrices

        :pdb: The pdb to process
        :yields: Yields a series of centers and rotation matrices
        for the corresponding tables
        """

        structure = self.structure(pdb)

        c = 0
        for residue in structure.residues():
            for name in residue.centers.definitions():
                center = residue.centers[name]
                if len(center) == 3:
                    yield mod.UnitCenters(unit_id=residue.unit_id(),
                                          name=name,
                                          pdb_id=pdb,
                                          x=float(center[0]),
                                          y=float(center[1]),
                                          z=float(center[2]))

            if hasattr(residue, 'rotation_matrix'):
                matrix = residue.rotation_matrix
                # if there are not enough atoms for the rotation matrix, it will be None
                if matrix is not None:
                    c += 1
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

