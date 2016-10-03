import operator as op
import itertools as it

from pymotifs import core
from pymotifs import models as mod

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

    def pdbs_in_nr(self, release, resolution='all', **kwargs):
        with self.session() as session:
            classes = mod.NrClasses
            chains = mod.NrChains
            ifes = mod.IfeInfo
            query = session.query(ifes.pdb_id).\
                join(chains, ifes.ife_id == chains.ife_id).\
                join(classes, classes.nr_class_id == chains.nr_class_id).\
                filter(chains.nr_release_id == release).\
                filter(classes.resolution == resolution).\
                distinct()
            return sorted({r.pdb_id for r in query})

    def to_process(self, pdbs, nr_release=None, **kwargs):
        if nr_release:
            return [tuple(self.pdbs_in_nr(nr_release, **kwargs))]
        if not pdbs:
            kwargs['all'] = True
        return super(Reporter, self).to_process(pdbs, **kwargs)

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
