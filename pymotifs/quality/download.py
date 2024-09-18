"""
Download and store validation reports. This will download all validation
reports from PDB and store them locally as a gzip file. In addition, it will
write an empty file for any report it could not download.
"""

import os
import ftplib
import urllib.request

import pymotifs.utils as ut
import pymotifs.core as core

import pymotifs.quality.utils as qut

from pymotifs.download import Writer


class Loader(core.Loader):
    """The loader to fetch and store quality data for structures.
    """
    dependencies = set()
    saver = Writer
    allow_no_data = True
    path = 'pub/pdb/validation_reports/{short}/{pdb}/{pdb}_validation.xml.gz'

    @property
    def ftp(self):
        """A FTP connection to the wwpdb.org FTP site.
        """
        if not hasattr(self, '_ftp'):
            self._ftp = ut.FTPFetchHelper('ftp.wwpdb.org')
        return self._ftp

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
        util = self._create(qut.Utils)
        return util.filename(pdb, compressed=True)

    def remote(self, pdb):
        """
        Get the path to validation report on PDB's FTP servers.

        Parameters
        ----------
        pdb : str
            The pdb id to use

        Returns
        -------
        path : str
            The path to the validation file.
        """
        return self.path.format(short=pdb[1:3].lower(), pdb=pdb.lower())

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

        Returns
        -------
        has_data : bool
            Check if there is data downloaded.
        """
        return os.path.exists(self.filename(pdb))

    def data(self, pdb, **kwargs):
        """
        Fetch the validation report from PDB and convert the entries there
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

        # Example: https://files.rcsb.org/pub/pdb/validation_reports/tn/4tna/4tna_validation.xml.gz
        pdb_url = "https://files.rcsb.org/pub/pdb/validation_reports/{short}/{pdb}/{pdb}_validation.xml.gz".format(short=pdb[1:3].lower(), pdb=pdb.lower())

        # Example: /usr/local/pipeline/hub-core/MotifAtlas/quality/validation-reports/4tna.xml.gz
        filename = self.filename(pdb)

        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.logger.info("Saving data in %s" % filename)

        try:
            # download from pdb_url and save the .gz file as filename
            urllib.request.urlretrieve(pdb_url, filename)
            return None
            # remote = self.remote(pdb)
            # return self.ftp(remote)
        except ftplib.error_perm:
            self.logger.warning("Could not fetch quality data for %s", pdb)
            return ''
        except ut.RetryFailedException:
            self.logger.warning("No quality found of %s", pdb)
            return ''
        except Exception as err:
            self.logger.exception(err)
            self.logger.warning("Unknown exception occurred")
            return ''
