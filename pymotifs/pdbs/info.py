"""Update pdb_info table. This uses PDB RESTful services for getting custom
reports. This will fetch the data for all the given PDBs if any one does not
have data. The old data will be overwritten as needed.
"""

from pymotifs import core
from pymotifs.utils.pdb import CustomReportHelper

from pymotifs import models as mod


class Loader(core.MassLoader):
    """This loads information about entire files into our database.
    """

    merge_data = True
    """We allow for merging data as the primary key is the pdb id."""

    dependencies = set()
    """This has no dependencies"""

    names = {
        'structureId': 'pdb_id',
        'structureTitle': 'title',
        'experimentalTechnique': 'experimental_technique',
        'depositionDate': 'deposition_date',
        'releaseDate': 'release_date',
        'revisionDate': 'revision_date',
        'ndbId': 'ndb_id',
        'resolution': 'resolution'
    }
    """This dict specifies the data to request from PDB's REST service (the
    keys) and what to call that data in our database (the values).
    """

    def rename(self, report):
        """Rename the columns of the report, and cover some values to the
        correct types. The columns will be renamed according to `Loader.names`
        and the resolution value will be converted to a float if it exists.

        Parameters
        ----------
        report : dict
            A dictonary with the keys in `self.names`.

        Returns
        -------
        renamed : dict
            The renamed dictonary
        """

        renamed = {}
        for key, name in self.names.items():
            renamed[name] = report.get(key)

        if renamed['resolution'] == '':
            renamed['resolution'] = None

        if renamed['resolution']:
            try:
                renamed['resolution'] = float(renamed['resolution'])
            except:
                renamed['resolution'] = None
                self.logger.error("Resoultion entry %s for %s is not a number",
                                  renamed['resolution'], renamed.get('pdb_id'))

        return renamed

    def has_data(self, pdb, **kwargs):
        """Check if we have data for the given PDB.

        Parameters
        ----------
        pdb : str
            The PDB to check

        Returns
        -------
        has_data : bool
            True if we have data for the PDB.
        """

        with self.session() as session:
            query = session.query(mod.PdbInfo).\
                filter_by(pdb_id=pdb).\
                limit(1)
            return bool(query.count())

    def data(self, pdbs, **kwargs):
        """Compute a report of the given PDBs. This will use PDB's custom REST
        service to compute a report for the given PDBs. This will then be
        processed slightly for storage in our database. The columns used in the
        report will be the keys in `Loader.names`.

        Parameters
        ----------
        pdbs : list
            A list of PDB ids

        Returns
        -------
        reports : list
            A list of PdbInfo objects to write to the database.
        """

        helper = CustomReportHelper(fields=self.names.keys())
        data = helper(pdbs)
        if not data:
            raise core.StageFailed("Could not load data for all pdbs %s" %
                                   str(pdbs))

        if len(data) != len(pdbs):
            self.logger.error("Could not get all data for all pdbs")

        return [mod.PdbInfo(**self.rename(report)) for report in data]
