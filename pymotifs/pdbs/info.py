"""
Update pdb_info and pdb_obsolete tables. Pdb_info only contains the current
pdb files, old files are removed. Uses PDB RESTful services for getting
custom reports.
"""

from pymotifs import core
from pymotifs.utils.pdb import CustomReportHelper

from pymotifs import models as mod


class Loader(core.MassLoader):
    """This loads information about entire files into our database.
    """

    merge_data = True
    dependencies = set()
    table = mod.PdbInfo

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

    def rename(self, report):
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
        with self.session() as session:
            query = session.query(mod.PdbInfo).\
                filter_by(pdb_id=pdb).\
                limit(1)
            return bool(query.count())

    def data(self, pdbs, **kwargs):
        helper = CustomReportHelper(fields=self.names.keys())
        data = helper(pdbs)
        if not data:
            raise core.StageFailed("Could not load data for all pdbs %s" %
                                   str(pdbs))

        if len(data) != len(pdbs):
            self.logger.error("Could not get all data for all pdbs")

        return [self.rename(report) for report in data]
