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

    # check for units that were missed and fill them in
    fill_in_missing = False

    def to_process(self, pdbs, **kwargs):
        """
        Do one query to determine which pdb ids do not have centers/rotations
        Since you need a center to have a rotation, just check for centers
        """

        if self.fill_in_missing:
            return pdbs

            # takes too much memory to try to find all unit ids that don't
            # have a center

        else:
            existing_ids = set()
            with self.session() as session:
                query = session.query(mod.UnitCenters.pdb_id).distinct()

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

        if self.fill_in_missing:
            # return an empty query so we process all pdb ids
            return session.query(mod.UnitCenters).filter_by(pdb_id='nonexistent_pdb_id')
        else:
            return session.query(mod.UnitCenters).filter_by(pdb_id=pdb)


    def data(self, pdb, **kwargs):
        """
        Load the structure, extract centers and rotation matrices

        :pdb: The pdb to process
        :yields: Yields a series of centers and rotation matrices
        for the corresponding tables
        """

        if self.fill_in_missing:
            # find unit ids in this pdb that have a center
            with self.session() as session:
                query = session.query(mod.UnitCenters).filter_by(pdb_id=pdb)
                has_center = set([u.unit_id for u in query])

            # find all unit ids in this pdb
            with self.session() as session:
                query = session.query(mod.UnitInfo.unit_id).filter_by(pdb_id=pdb)
                unit_ids = set([u.unit_id for u in query])

            if len(unit_ids - has_center) == 0:
                raise core.Skip("All unit ids in %s already have centers" % pdb)

        else:
            has_center = set()

        structure = self.structure(pdb)

        c = 0
        for residue in structure.residues():
            # if the unit has any center at all, skip it
            # if we ever have to add some new center, this logic will need to change

            if residue.unit_id() in has_center:
                continue

            for name in residue.centers.definitions():
                center = residue.centers[name]

                if self.fill_in_missing:
                    self.logger.info("Adding center for %s value %s" % (residue.unit_id(), center))
                    print("Adding center for %s value %s" % (residue.unit_id(), center))

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

