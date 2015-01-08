import os

import core
from .cif import CifStructure

from fr3d.cif.writer.cifatom import CifAtom


class Exporter(core.PdbExporter):
    paramters = (CifStructure)

    def filename(self, pdb):
        return os.path.join(self.config['fr3d_root'], "PDBFiles",
                            pdb + ".cifatoms")

    def process(self, structure):
        with open(self.filename(structure.pdb), 'rb') as out:
            writer = CifAtom(out)
            writer(structure)
