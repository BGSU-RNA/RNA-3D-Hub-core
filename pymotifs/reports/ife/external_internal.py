import operator as op
import itertools as it

from pymotifs import core

from pymotifs.ife.helpers import IfeLoader


class Reporter(core.Reporter):
    headers = [
        'Pdb',
        'Chain1',
        'Chain2',
        'Length',
        'Internal',
        'External',
    ]

    def ifes(self, pdb):
        loader = IfeLoader(self.config, self.session)
        try:
            ifes, interactions = loader(pdb)
        except Exception:
            return []
        data = []

        structured = it.ifilter(op.attrgetter('is_structured'), ifes)
        for ife1, ife2 in it.combinations(structured, r=2):
            data.append({
                'Pdb': pdb,
                'Chain1': ife1.chain,
                'Chain2': ife2.chain,
                'Length': min(ife1.length, ife2.length),
                'Internal': min(ife1.internal, ife2.internal),
                'External': interactions[ife1.chain][ife2.chain],
            })
        return data

    def data(self, pdbs, **kwargs):
        data = []
        for pdb in pdbs:
            data.extend(self.ifes(pdb))
        return data
