"""This is a module to export a cif file into a cifatoms format for matlab. It
writes the cifatoms files to the
"""

import os
from contextlib import contextmanager

from pymotifs import core
from pymotifs import savers

from pymotifs.download import Downloader

from fr3d.cif.writer import CifAtom
from fr3d.cif.reader import ComplexOperatorException


class CifAtomExportFailed(Exception):
    """Raised when the CifAtom export fails. There are some cases when this
    fails silently while writing out a file. We want to raise an exception so
    we know things fail.
    """
    pass


class CifAtomSaver(savers.FileHandleSaver):
    """A Saver to write Cifatom files. This will write a Structure into a
    cifatom file.
    """

    allow_empty = False

    @contextmanager
    def writer(self, pdb, **kwargs):
        with super(self, CifAtomSaver).writer(pdb, **kwargs) as handle:
            writer = CifAtom(handle)
            try:
                yield writer
            except ComplexOperatorException:
                self.logger.warning("Cannot export %s to cifatom", pdb)


class Exporter(core.Loader):
    """Will export files from the mmCIF format to a cifatom format readable by
    matlab programs.
    """

    dependencies = set([Downloader])
    saver = CifAtomSaver

    def filename(self, pdb, **kwargs):
        return os.path.join(self.config['locations']['fr3d_root'], "PDBFiles",
                            pdb + ".cifatoms")

    def is_missing(self, entry, **kwargs):
        """Will check if the file produce by filename() exists.
        """
        filename = self.filename(entry)
        try:
            return not os.path.isfile(filename) and \
                os.path.getsize(filename) > 0
        except:
            return False

    def data(self, pdb, **kwargs):
        """Will load the structure for the given PDB id.
        """
        return self.structure(pdb)
