"""A loader to store the coordinates for each residue in the database. This
will write one row per non water residue into the unit_coordinates table in the
database. Data is written in the CIF format. This will write only the ATOM
level entries (atom_site in cif) but will not include the header lines like
'loop_' or '_atom_site.group_PDB'.
"""

from cStringIO import StringIO

import pymotifs.core as core
from pymotifs import models as mod

from fr3d.cif.writer import CifAtom
from fr3d.data import Structure

from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    """The loader to store unit_coordinates data.
    """

    dependencies = set([InfoLoader])

    def query(self, session, pdb):
        """Create a query to find all all entries in `units_coordinates` for
        the given PDB id.

        Parameters
        ----------
        session : pymotifs.core.Session
            The session to use.

        pdb : str
            The PDB id to use.

        Returns
        -------
        query : Query
            The query for the given structure
        """

        return session.query(mod.UnitCoordinates).\
            join(mod.UnitInfo,
                 mod.UnitInfo.unit_id == mod.UnitCoordinates.unit_id).\
            filter(mod.UnitInfo.pdb_id == pdb)

    def coordinates(self, pdb, residue):
        """Compute a string of the coordinates in CIF format (the atom_site
        block) for the given residue. This will exclude the header and trailing
        lines that are part of the atom_site entries. This is because these are
        entries are meant to be concatenated together for the coordinate server
        later.

        Parameters
        ----------
        pdb : str
            The PDB id to use

        residue : fr3d.data.Component
            The residue to convert

        Returns
        -------
        coordinates : str
            A string that represents CIF data formated data for the given
            residue.
        """

        structure = Structure([residue], pdb=pdb)
        sio = StringIO()
        writer = CifAtom(sio, unit_ids=False)
        writer(structure)
        raw = sio.getvalue()
        # Lines that star with data_, loop_ and _ or # are header/comment lines
        # and so are excluded
        coords = []
        for line in raw.split('\n'):
            if not line or \
                    line.startswith('data_') or \
                    line.startswith('loop_') or \
                    line[0] in set('_#'):
                continue
            coords.append(line)
        return '\n'.join(coords)

    def data(self, pdb, **kwargs):
        """Compute the coordinate entries for the given PDB. This will exclude
        water molecules as those aren't generally worth displaying in
        coordinate server.

        Parameters
        ----------
        pdb : str
            The PDB id to use

        Yields
        ------
        coord : UnitCoordinates
            A UnitCoordinates object with the coordinates to write.
        """

        structure = self.structure(pdb)
        for unit in structure.residues():
            if unit.sequence == 'HOH':
                continue
            coord = self.coordinates(pdb, unit)
            if not coord:
                raise core.InvalidState("No coordinates computed for %s" %
                                        unit)

            yield mod.UnitCoordinates(
                unit_id=unit.unit_id(),
                coordinates=coord,
            )
