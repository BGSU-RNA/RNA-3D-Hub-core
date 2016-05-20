"""A loader to store the coordinates for each resideu in the database. This
will write one row per non water residue into the unit_coordinates table in the
database. Data is written in the CIF format. This will write only the ATOM
level entries (atom_site in cif) but will not includ the header lines like
'loop_' or '_atom_site.group_PDB'.
"""

from cStringIO import StringIO

import pymotifs.core as core

from fr3d.cif.writer import CifAtom
from fr3d.data import Structure

from pymotifs.models import UnitInfo
from pymotifs.models import UnitCoordinates
from pymotifs.units.info import Loader as InfoLoader


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader])
    table = UnitCoordinates

    def query(self, session, pdb):
        return session.query(UnitCoordinates).\
            join(UnitInfo, UnitInfo.unit_id == UnitCoordinates.unit_id).\
            filter(UnitInfo.pdb_id == pdb)

    def coordinates(self, pdb, residue):
        structure = Structure([residue], pdb=pdb)
        sio = StringIO()
        writer = CifAtom(sio, units=False)
        writer(structure)
        raw = sio.getvalue()
        # Skip the first 24 lines which consist of all the header information
        # for cif atoms.
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
        structure = self.structure(pdb)
        for unit in structure.residues():
            if unit.sequence == 'HOH':
                continue
            yield {
                'unit_id': unit.unit_id(),
                'coordinates': self.coordinates(pdb, unit)
            }
