"""
Update pdb_info and pdb_obsolete tables. Pdb_info only contains the current
pdb files, old files are removed. Uses PDB RESTful services for getting
custom reports.
"""

import core
from utils.pdb import CustomReportHelper

from models import PdbInfo


class PdbInfoLoader(core.Loader):
    """This loads information about entire files into our database.
    """

    fields = [
        'structureId',
        'structureTitle',
        'experimentalTechnique',
        'depositionDate',
        'releaseDate',
        'revisionDate',
        'ndbId',
        'resolution',
        'classification',
    ]

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
        super(PdbInfoLoader, self).__init__(*args)
        self.helper = CustomReportHelper(fields=self.fields)

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(PdbInfo).filter_by(id=pdb)
            return bool(query.count())

    def remove(self, pdb):
        with self.session() as session:
            session.query(PdbInfo).\
                filter_by(id=pdb).\
                delete(synchronize_session='fetch')

    def rename(self, report):
        renamed = {}
        for key, name in self.names.items():
            renamed[name] = report[key]
        return renamed

    def data(self, pdb):
        return [PdbInfo(**self.rename(report)) for report in self.helper(pdb)]
