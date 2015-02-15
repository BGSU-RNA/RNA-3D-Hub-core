from Bio.Alphabet import ThreeLetterProtein

from pymotifs import core
from pymotifs.models import PdbUnitIdCorrespondence

AA = [seq.upper() for seq in ThreeLetterProtein().letters]


class TranslationFailed(Exception):
    """This is raised when we could not translate an old style id to a new
    style id.
    """


def component_type(unit):
    seq = unit.sequence.upper()
    if seq in ['A', 'C', 'G', 'U']:
        return 'rna'
    if seq == 'HOH':
        return 'water'
    if seq in AA:
        return 'aa'
    if seq in ['DA', 'DC', 'DG', 'DT']:
        return 'dna'
    return None


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
