from Bio.Alphabet import ThreeLetterProtein

from pymotifs import core
from pymotifs import models as mod

AA = [seq.upper() for seq in ThreeLetterProtein().letters]


class TranslationFailed(Exception):
    """This is raised when we could not translate an old style id to a new
    style id.
    """
    pass


def component_type(unit):
    seq = unit.sequence.upper()
    if seq in ['A', 'C', 'G', 'U']:
        return 'rna'
    if seq in ['DA', 'DC', 'DG', 'DT']:
        return 'dna'
    if seq == 'HOH':
        return 'water'
    if seq in AA:
        return 'aa'
    ions = unit.unit_id().split('|')[3].upper()
    if ions in ['CU', 'FE', 'MG', 'NI', 'MN', 'K', 'NA', 'MO', 'CO', 'ZN', 'W', 'CA', 'V']:
        return 'ion'
    nt = unit.type.lower()
    if nt == 'rna linking':
        return 'rna'
    if nt == 'dna linking':
        return 'dna'
    if nt == 'l-peptide linking':
        return 'aa'
    return None


class Translator(object):
    def __init__(self, session):
        self.session = core.Session(session)

    def translate(self, raw):
        mapping = {}
        corr = mod._PdbUnitIdCorrespondence
        with self.session() as session:
            query = session.query(corr.unit_id, corr.old_id)

            nt_ids = raw
            if isinstance(raw, str):
                nt_ids = raw.split(',')

            query = query.filter(corr.old_id.in_(nt_ids))

            for result in query:
                mapping[result.old_id] = result.unit_id

        units = []
        for nt_id in nt_ids:
            if nt_id not in mapping:
                raise TranslationFailed("Could not find %s" % nt_id)
            units.append(mapping[nt_id])

        return units
