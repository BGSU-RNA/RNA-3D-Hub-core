"""
Update pdb_info and pdb_obsolete tables. Pdb_info only contains the current
pdb files, old files are removed. Uses PDB RESTful services for getting
custom reports.
"""

from pymotifs import core
from pymotifs.utils.pdb import CustomReportHelper

from pymotifs.models import PdbInfo


class Loader(core.MassLoader):
    """This loads information about entire files into our database.
    """
    merge_data = True

    names = {
        'structureId': 'id',
        'structureTitle': 'title',
        'experimentalTechnique': 'experimental_techinque',
        'depositionDate': 'deposition_date',
        'releaseDate': 'release_date',
        'revisionDate': 'revision_date',
        'ndbId': 'ndb_id',
        'resolution': 'resolution'
    }

    def __init__(self, *args):
        super(Loader, self).__init__(*args)
        self.helper = CustomReportHelper(fields=self.names.keys())

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
                                  renamed['resolution'], renamed['id'])

        return renamed

    def data(self, pdbs, **kwargs):
        data = self.helper(pdbs)
        if not data:
            raise core.StageFailed("Could not load data for all pdbs %s" %
                                   pdbs)
        return [PdbInfo(**self.rename(report)) for report in data]
