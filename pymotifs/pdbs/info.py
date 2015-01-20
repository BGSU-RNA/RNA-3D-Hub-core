"""
Update pdb_info and pdb_obsolete tables. Pdb_info only contains the current
pdb files, old files are removed. Uses PDB RESTful services for getting
custom reports.
"""

from pymotifs import core
from pymotifs.utils.pdb import CustomReportHelper

from pymotifs.models import PdbInfo


class Loader(core.SimpleLoader):
    """This loads information about entire files into our database.
    """

    names = {
        'structureId': 'id',
        'structureTitle': 'title',
        'experimentalTechnique': 'experimental_technique',
        'depositionDate': 'deposition_date',
        'releaseDate': 'release_date',
        'revisionDate': 'revision_date',
        'ndbId': 'ndb_id',
        'resolution': 'resolution',
        'classification': 'classification',
    }

    def __init__(self, *args):
        super(Loader, self).__init__(*args)
        self.helper = CustomReportHelper(fields=self.names.keys())

    def query(self, session, pdb):
        return session.query(PdbInfo).filter_by(id=pdb)

    def rename(self, report):
        renamed = {}
        for key, name in self.names.items():
            renamed[name] = report[key]
        return renamed

    def data(self, pdb):
        return [PdbInfo(**self.rename(report)) for report in self.helper(pdb)]