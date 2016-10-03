from pymotifs import core

from pymotifs.ife.helpers import IfeLoader


class Reporter(core.Reporter):
    headers = [
        'Pdb',
        'Chain1',
        'Chain2',
        'Internal',
        'External',
    ]

    def interactions(self, pdb):
        loader = IfeLoader(self.config, self.session)
        try:
            ifes, interactions = loader(pdb)
        except Exception:
            return []
        data = []
        for ife1 in ifes:
            if not ife1.is_structured:
                continue
            for ife2 in ifes:
                data.append({
                    'Pdb': pdb,
                    'Chain1': ife1.chain,
                    'Chain2': ife2.chain,
                    'Internal': min(ife1.internal, ife2.internal),
                    'External': interactions[ife1.chain][ife2.chain],
                })
        return data

    def data(self, pdbs, **kwargs):
        data = []
        for pdb in pdbs:
            data.extend(self.interactions(pdb))
        return data
