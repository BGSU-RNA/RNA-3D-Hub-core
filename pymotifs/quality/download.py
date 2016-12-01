"""Download and store validation reports. This will download all validation
reports from PDB and store them locally as a gzip file. In addition, it will
write an empty file for any report it could not download.
"""

import os

import pymotifs.utils as ut
import pymotifs.core as core

import pymotifs.quality.utils as qut


class Loader(core.Loader):
    """The loader to fetch and store quality data for structures.
    """

    dependencies = set()
    saver = core.FileHandleSaver
    path = 'pub/pdb/validation_reports/{short}/{pdb}/{pdb}_validation.xml.gz'

    def filename(self, pdb, **kwargs):
        """Get the filename where data for this PDB should be stored.

        Parameters
        ----------
        pdb : str
            The PDB id to use.

        Returns
        -------
        path : str
            The path to write to.
        """
        return qut.filename(self.config, pdb, compressed=True)

    def remove(self, entry, **kwargs):
        """Remove the data if it exists.

        Parameters
        ----------
        pdb : str
            The PDB to remove data for.
        """
        if self.has_data(entry) and not kwargs.get('dry_run'):
            os.remove(self.filename(entry))

    def has_data(self, pdb, **kwargs):
        """Check if this has data. This is done by checking if stored file
        exists.

        Parameters
        ----------
        pdb : str
            The pdb to lookup.
        ignore_empty : bool, False
            This will ingore

        Returns
        -------
        has_data : bool
            Check if there is data downloaded.
        """
        return os.path.exists(self.filename(pdb))

    def data(self, pdb, **kwargs):
        """Compute the quality assignments for residues in the structure. This
        will fetch the validation report from PDB and convert the entries there
        into forms suitable to write to the database. If the report has no RSR
        or DCC data then a `core.Skip` exception will be raised.

        Parameters
        ----------
        pdb : str
            The pdb id to use.

        Returns
        -------
        data: iterable
            An iterable of a quality assignments to store in the database.
        """

        filename = self.filename(pdb)
        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)

        try:
            fetcher = ut.FTPFetchHelper('ftp.wwpdb.org')
            return fetcher(filename)
        except ut.RetryFailedException:
            self.logger.warning("No quality found of %s", pdb)
            return ''
