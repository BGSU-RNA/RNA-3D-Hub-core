from pymotifs import core
from pymotifs.models import PdbUnitIdCorrespondence


class TranslationFailed(Exception):
    """This is raised when we could not translate an old style id to a new
    style id.
    """


class Translator(object):
    def __init__(self, session):
        self.session = core.Session(session)

    def translate(self, raw):
        mapping = {}
        with self.session() as session:
            query = session.query(PdbUnitIdCorrespondence.unit_id,
                                  PdbUnitIdCorrespondence.old_id)

            nt_ids = raw
            if isinstance(raw, str):
                nt_ids = raw.split(',')

            query = query.filter(PdbUnitIdCorrespondence.old_id.in_(nt_ids))

            for result in query:
                mapping[result.old_id] = result.unit_id

        units = []
        for nt_id in nt_ids:
            if nt_id not in mapping:
                raise TranslationFailed("Could not find %s" % nt_id)
            units.append(mapping[nt_id])

        return units
