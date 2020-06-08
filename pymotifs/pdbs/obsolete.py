"""A stage to load all obsolete data from PDB. This will use PDB's FTP site to
fetch the current list of obsolete PDB ids and store them in our database. Note
that running this on partially setup database (like the one produced by
bootstrapping) will cause foreign key issues. This is because we get all PDB's
that have ever been obsoleted and a partially setup database will not have data
on all PDB's.
"""

import time
from datetime import datetime
from cStringIO import StringIO

from pymotifs import core
from pymotifs import utils
from pymotifs import models as mod

from pymotifs.pdbs.info import Loader as InfoLoader

import urllib2


class Parser(object):
    """A class to help parsing the file fetched from PDB's FTP site. This is a
    callable object that will turn the large chunk of text into a list of
    dictonaries.
    """

    def __call__(self, text):
        """Parse the given text.

        Parameters
        ----------
        text : str
            The text to parse

        Returns
        -------
        obsolete : list
            A list of dictonaries that contains a 'pdb_obsolete_id', 'date',
            and 'replaced_by' keys.
        """

        data = []
        for line in StringIO(text).readlines():
            # OBSLTE    26-SEP-06 2H33     2JM5 2OWI
            if 'OBSLTE' in line.split():
                parts = line.split()
                obsolete_date = datetime.strptime(parts[1], '%d-%b-%y')
                replaced = ",".join(parts[3:])
                if not replaced:
                    replaced = None
                data.append({
                    'pdb_obsolete_id': parts[2],
                    'date': obsolete_date,
                    'replaced_by': replaced
                })
        return data


class Loader(core.MassLoader):
    """The actual loader.
    """

    merge_data = True
    """We can just update data in place, so we do."""

    dependencies = set([InfoLoader])
    """We only depend on the pdbs.info loader."""

    @property
    def table(self):
        return mod.PdbObsolete

    def has_data(self, pdb_id, **kwargs):
        """We always return False so we always fetch new updates.

        Parameters
        ----------
        pdb_id : str
            Ignored.

        Returns
        -------
        False
        """
        return False

    def data(self, *args, **kwargs):
        """Download the file with all obsolete structures over ftp and parse
        them into a storable format.

        Returns
        -------
        data : list
            A list of dictonaries as from `Parser`.
        """
        attempts = 0

        while attempts < 100:
            try:
                response = urllib2.urlopen('https://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat')
                html = response.read()
                return html
            except Exception as err:
                attempts += 1
                print("Failed %d times to get obsolete ids via URL" % attempts)
                time.sleep(5)

        attempts = 0

        while attempts < 100:
            try:
                ftp = utils.FTPFetchHelper('ftp.wwpdb.org', parser=Parser())
                return ftp('/pub/pdb/data/status/obsolete.dat')
            except Exception as err:
                attempts += 1
                print("Failed %d times to get obsolete ids" % attempts)
                time.sleep(5)
                self.logger.critical("Could not get obsolete ids")
                self.logger.exception(err)

        raise core.StageFailed("Could not get obsolete ids")
