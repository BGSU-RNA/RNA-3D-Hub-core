"""This is a module to export a cif file into a cifatoms format for matlab. It
writes the cifatoms files to the
"""

import os
from contextlib import contextmanager

from pymotifs import core

from pymotifs.download import Downloader

from fr3d.cif.writer import CifAtom
from fr3d.cif.reader import ComplexOperatorException


class CifAtomExportFailed(Exception):
    """Raised when the CifAtom export fails. There are some cases when this
    fails silently while writing out a file. We want to raise an exception so
    we know things fail.
    """
    pass


class CifAtomSaver(core.FileHandleSaver):
    """A Saver to write Cifatom files. This will write a Structure into a
    cifatom file.
    """

    allow_empty = False

    @contextmanager
    def writer(self, pdb, **kwargs):
        with super(CifAtomSaver, self).writer(pdb, **kwargs) as handle:
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
        """Create the filename for the given pdb file.
        """
        return os.path.join(self.config['locations']['fr3d_root'], "PDBFiles",
                            pdb + ".cifatoms")

    def has_data(self, entry, **kwargs):
        """Will check if the file produce by filename() exists and is not
        empty.
        """
        filename = self.filename(entry)
        try:
            return os.path.isfile(filename) and \
                os.path.getsize(filename) > 0
        except:
            return False

    def remove(self, pdb, **kwargs):
        """Remove the file for the given pdb.
        """
        if kwargs.get('dry_run'):
            self.logger.info("Not doing anything, is dry run")
            return None

        filename = self.filename(pdb)
        if os.path.isfile(filename):
            os.remove(filename)

    def data(self, pdb, **kwargs):
        """Will load the structure for the given PDB id.
        """
        return self.structure(pdb)
