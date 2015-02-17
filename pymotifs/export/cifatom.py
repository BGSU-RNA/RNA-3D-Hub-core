"""This is a module to export a cif file into a cifatoms format for matlab. It
writes the cifatoms files to the
"""

import os

from pymotifs import core

from fr3d.cif.writer import CifAtom


class Exporter(core.Stage):

    def filename(self, pdb):
        return os.path.join(self.config['locations']['fr3d_root'], "PDBFiles",
                            pdb + ".cifatoms")

    def process(self, pdb, **kwargs):
        with open(self.filename(pdb), 'wb') as out:
            writer = CifAtom(out)
            writer(self.structure(pdb))
