"""Parse and store validation reports for pdb level information. This will
process the downloaded validation reports and use them to populate the
units_quality table.
"""

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.quality.utils import Parser
from pymotifs.quality.utils import known

from pymotifs.pdb.info import Loader as PdbLoader
from pymotifs.quality.download import Loader as Downloader


class Loader(core.SimpleLoader):
    """The loader to fetch and store quality data for structures.
    """

    dependencies = set([PdbLoader, Downloader])

    mapping = {
        'percent-RSRZ-outliers': float,
        'relative-percentile-percent-RSRZ-outliers': float,
    }

    def to_process(self, pdbs, **kwrags):
        return sorted(set(pdbs) - set(known(has_data=False)))

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
        data : iterable
            An iterable of a quality assignments to store in the database.
        """

        with open(self.filename(pdb), 'rb') as raw:
            parser = Parser(raw.read())
        entity = parser.entity()

        data = {'pdb_id': pdb}
        for key, fn in self.mapping.items():
            value = entity.get(key, None)
            if value is not None:
                value = fn(value)
            data[key.replace('-', '_')] = value
        return mod.PdbQuality(**data)
