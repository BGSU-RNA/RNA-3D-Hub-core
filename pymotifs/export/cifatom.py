"""This is a module to export a cif file into a cifatoms format for matlab. It
writes the cifatoms files to the
"""

import os

from pymotifs import core

from pymotifs.download import Downloader

from fr3d.cif.writer import CifAtom
from fr3d.cif.reader import ComplexOperatorException


class Exporter(core.Stage):
    dependencies = set([Downloader])

    def filename(self, pdb):
        return os.path.join(self.config['locations']['fr3d_root'], "PDBFiles",
                            pdb + ".cifatoms")

    def is_missing(self, entry, **kwargs):
        """Will check if the file produce by filename() exists.
        """
        return not os.path.exists(self.filename(entry))

    def process(self, pdb, **kwargs):
        with open(self.filename(pdb), 'wb') as out:
            writer = CifAtom(out)
            try:
                writer(self.structure(pdb))
            except ComplexOperatorException:
                self.logger.warning("Cannot export %s to cifatom", pdb)
