
from pymotifs import core
from pymotifs import models as mod

class TranslationFailed(Exception):
    """
    This is raised when we could not translate an old style id to a new
    style id.
    """
    pass


def component_type(residue,unit=None):
    """
    This is the most complete way to determine the unit_type_id.
    Put all your effort here.
    residue is a full unit as read from a .cif file
    unit is the name of the unit, just a sequence
    """

    if residue:
        seq = residue.sequence.upper()
    else:
        seq = unit.upper()

    if seq in ['A','C','G','U']:
        return 'rna'
    if seq in ['DA','DC','DG','DT']:
        return 'dna'
    if seq == 'HOH':
        return 'water'

    # most common amino acids
    AA = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","SEC","PYL"]
    if seq in AA:
        return 'aa'

    ions = set(['AG','AU','BA','BR','CA','CD','CL','CO','CR','CS','CU','F','FE','GD','HG','IR','K','LI','LU','MG','MN','MO','NA','NI','O','OH','OS','PB','PT','RB','RE','RU','SM','SR','TB','TL','W','V','YB','ZN'])
    if seq in ions:
        return 'ion'

    if residue:
        unit_type = residue.type.lower()
        if unit_type == 'rna linking':
            return 'rna'
        if unit_type == 'rna oh 3 prime terminus':
            return 'rna'
        if unit_type == 'l-rna linking':
            return 'rna'
        if unit_type == 'dna linking':
            return 'dna'
        if unit_type == 'dna oh 3 prime terminus':
            return 'dna'
        if unit_type == 'l-dna linking':
            return 'dna'
        if unit_type == 'l-peptide linking':
            return 'aa'

    from fr3d.modified.mapping import modified_base_to_parent
    if seq in modified_base_to_parent:
        if modified_base_to_parent[seq] in ['A','C','G','U']:
            return 'rna'
        else:
            return 'dna'

    return None


class Translator(object):
    def __init__(self,session):
        self.session = core.Session(session)

    def translate(self,raw):
        mapping = {}
        corr = mod._PdbUnitIdCorrespondence
        with self.session() as session:
            query = session.query(corr.unit_id,corr.old_id)

            nt_ids = raw
            if isinstance(raw,str):
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
