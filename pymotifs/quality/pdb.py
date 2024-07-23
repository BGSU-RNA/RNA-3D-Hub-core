"""Parse and store validation reports for pdb level information. This will
process the downloaded validation reports and use them to populate the
pdb_quality table.
"""
import os

import os

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.quality.utils import Parser
from pymotifs.quality.utils import Utils

from pymotifs.pdbs.info import Loader as PdbLoader
from pymotifs.quality.download import Loader as Downloader


class Loader(core.SimpleLoader):
    """The loader to fetch and store quality data for structures.
    """

    dependencies = set([PdbLoader, Downloader])

    def to_process(self, pdbs, **kwargs):
        """Get the pdbs to process. These will be the structures that have an
        non-empty validation report downloaded.

        Parameters
        ----------
        pdbs : list
            The pdbs to get data for

        Returns
        -------
        pdbs : list
            The list of PDBs to process.
        """
        util = self._create(Utils)
        return sorted(set(pdbs) - set(util.known(has_data=False)))

    def query(self, session, pdb):
        """Generate a query to find all entries in units_quality for the given
        PDB id.

        Attributes
        ----------
        session : Session
            The `Session` to use.

        pdb : str
            The PDB id to use.

        Returns
        -------
        query : Query
            Returns an SqlAlchemy query for all entires in units_quality with
            the given PDB id.
        """
        return session.query(mod.PdbQuality).\
            filter_by(pdb_id=pdb)

    def parse(self, filename):
        """Parse the file to extract the structure level data.

        Parameters
        ----------
        filename : str
            The file to parse.

        Returns
        -------
        data : mod.UnitQuality
            The quality data for the structure.
        """
        if not os.path.exists(filename):
            raise core.Skip("Missing file %s" % filename)
        with open(filename, 'rb') as raw:
            parser = Parser(raw.read())
        entity = parser.entity()
        return mod.PdbQuality(**entity)

    def data(self, pdb, **kwargs):
        """Compute the quality assignments for the structure.

        Parameters
        ----------
        pdb : str
            The pdb id to use.

        Returns
        -------
        data : mod.UnitQuality
            The quality data for the structure.
        """
        filename = self._create(Utils).filename(pdb)
        if not os.path.exists(filename):
            raise core.Skip("No quality for %s" % pdb)
        return self.parse(filename)
