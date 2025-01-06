"""
A loader to store the coordinates for each residue in the database.
This will write one row per atom into the unit_coordinates table in the
database. Data is written in the CIF format. This will write only the ATOM
level entries (atom_site in cif) but will not include the header lines like
'loop_' or '_atom_site.group_PDB'.

Ligands, ions, and water molecules were not being written until January 2025.

"""

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.units.info import Loader as InfoLoader

from fr3d.cif.writer import CifAtom
from fr3d.data import Structure

import sys

# safe for python 2 and 3
if sys.version_info[0] == 2:
    from cStringIO import StringIO
else:
    from io import StringIO

class Loader(core.SimpleLoader):
    """
    The loader to store unit_coordinates data.
    """

    dependencies = set([InfoLoader])

    allow_no_data = True

    # check for units that were missed and fill them in
    fill_in_missing = False

    def query(self, session, pdb):
        """
        Create a query to find all entries in `units_coordinates` for the
        given PDB id.

        Parameters
        ----------
        session : pymotifs.core.Session
            The session to use.

        pdb : str
            The PDB id to use.

        Returns
        -------
        query : Query
            The query for the given structure.
        """

        if self.fill_in_missing:
            # return an empty query so we process all pdb ids
            return session.query(mod.UnitCoordinates).filter_by(unit_id='nonexistent_unit_id')
        else:
            return session.query(mod.UnitCoordinates).\
                join(mod.UnitInfo,
                    mod.UnitInfo.unit_id == mod.UnitCoordinates.unit_id).\
                filter(mod.UnitInfo.pdb_id == pdb)


    def coordinates(self, pdb, residue):
        """
        Compute a string of the coordinates in CIF format (the atom_site
        block) for the given residue. Exclude the header and trailing lines
        that are part of the atom_site entries, because these entries are meant
        to be concatenated together for the coordinate server later.

        Parameters
        ----------
        pdb : str
            The PDB id to use.

        residue : fr3d.data.Component
            The residue to convert.

        Returns
        -------
        coordinates : str
            A string that represents CIF-formatted data for the given residue.
        """

        # make the given residue into a structure
        structure = Structure([residue], pdb=pdb)
        sio = StringIO()
        writer = CifAtom(sio, unit_ids=False, protect_lists_of_lists=True)
        writer(structure)
        raw = sio.getvalue()
        coords = []
        for line in raw.split('\n'):
            # Exclude header/comment lines that start with: 1) "data_",
            # 2) "loop_",  3) "_", or 4) "#".
            if not line or \
                    line.startswith('data_') or \
                    line.startswith('loop_') or \
                    line[0] in set('_#'):
                continue
            line = line.replace("    "," ").replace("  "," ").replace("  "," ")
            coords.append(line)
        return '\n'.join(coords)


    def data(self, pdb, **kwargs):
        """
        Compute the coordinate entries for the given PDB.
        This will exclude water molecules as those aren't generally
        worth displaying in the coordinate server.

        Parameters
        ----------
        pdb : str
            The PDB id to use.

        Yields
        ------
        coord : UnitCoordinates
            A UnitCoordinates object with the coordinates to write.
        """

        if self.fill_in_missing:
            # find all unit ids in this pdb
            with self.session() as session:
                query = session.query(mod.UnitInfo.unit_id).filter_by(pdb_id=pdb)
                unit_ids = set([u.unit_id for u in query])

            # find unit ids in this pdb that have coordinates
            with self.session() as session:
                query = session.query(mod.UnitCoordinates).\
                                join(mod.UnitInfo,
                                    mod.UnitInfo.unit_id == mod.UnitCoordinates.unit_id).\
                                filter(mod.UnitInfo.pdb_id == pdb)
                has_coordinates = set([u.unit_id for u in query])

            if len(unit_ids - has_coordinates) == 0:
                raise core.Skip("All unit ids in %s already have centers" % pdb)

        else:
            has_coordinates = set()


        # read the .cif file
        structure = self.structure(pdb)

        # loop over units in the structure, RNA, DNA, protein, ligands, ions, etc.
        for unit in structure.residues():
            if unit.unit_id() in has_coordinates:
                continue

            # pass the pdb id and the entire unit to self.coordinates
            coord = self.coordinates(pdb, unit)
            self.logger.debug("data: PDB: %s" % pdb)
            self.logger.debug("data: unit: %s" % unit)
            self.logger.debug("data: coordinates: %s" % coord)
            if not coord:
                raise core.InvalidState("No coordinates computed for %s" % unit)

            if self.fill_in_missing:
                print("Adding coordinates for %s value %s" % (unit.unit_id(), coord))

            yield mod.UnitCoordinates(
                unit_id=unit.unit_id(),
                coordinates=coord,
            )
