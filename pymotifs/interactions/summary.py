from pymotifs import core

from pymotifs.models import PdbBestChainsAndModels
from pymotifs.models import PdbHelixLoopInteractionSummary as Summary
from pymotifs.interactions.pairwise import Loader as InterLoader

QUERY_TEMPLATE = '''
SELECT
    count(if(f_lwbp regexp '^[ct]' {unique}, 1, NULL)) as 'bps',
    count(if(f_stacks regexp '^s' {unique}, 1, NULL)) as 'stacks',
    count(if(f_bphs regexp '^[1-9]'
            or (f_bphs = '0BPh' and iPdbSig != jPdbSig), 1, NULL)) as 'bphs'
FROM pdb_pairwise_interactions
JOIN {table1} as N1
ON
    iPdbSig = N1.id
JOIN {table2} as N2
ON
    jPdbSig = N2.id
WHERE
    pdb_id = :pdb
    AND N1.chain = N2.chain
    AND N1.chain = :chain
    AND (
        (f_lwbp regexp '^[ct]' {unique}) or
        (f_stacks regexp '^s' {unique}) or
        (f_bphs regexp '^[1-9]') or
        (f_bphs = '0BPh' and iPdbSig != jPdbSig)
    )
    AND f_crossing {operator} 4
;
'''


class Loader(core.Loader):
    name = 'interaction_summary'
    update_gap = False
    dependencies = set(InterLoader)

    def has_data(self, entry):
        pdb, chain = entry
        with self.session() as session:
            query = session.query(Summary).filter_by(pdb=pdb, chain=chain)
            return bool(query.count())

    def transform(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(PdbBestChainsAndModels).filter_by(pdb_id=pdb)

            best = query.first()
            if not best:
                raise core.SkipPdb("Could not get best chain for %s" % pdb)

            return [(pdb, chain) for chain in best.best_chains.split(',')]

    def remove(self, entry):
        pdb, chain = entry
        with self.session() as session:
            session.query(Summary).filter_by(pdb=pdb, chain=chain).\
                delete(synchronize_session='fetch')

    def table(self, element):
        if element == 'helix':
            return 'nt_helix_view'
        if element == 'loop':
            return 'nt_loop_view'
        raise ValueError("Unknown element type %s" % element)

    def build(self, element1, element2, range_type):
        if range_type == 'sr':
            operator = '<'
        elif range_type == 'lr':
            operator = '>='
        else:
            raise ValueError("Unknown range type: %s" % range_type)

        unique = ''
        if element1 == element2:
            unique = 'and iPdbSig < jPdbSig'

        return QUERY_TEMPLATE.format(unique=unique, operator=operator,
                                     table1=self.table(element1),
                                     table2=self.table(element2))

    def query(self, pdb, chain, element1, element2, range_type):
        query = self.build(element1, element2, range_type)
        with self.session() as session:
            results = session.execute(query, {'pdb': pdb, 'chain': chain})

            results = results.fetchone()
            if not results:
                raise core.SkipValue("Couldn't compute %s %s" %
                                     element1, element2)

            return results['bps'], results['stacks'], results['bphs']

    def data(self, entry, **kwargs):
        pdb, chain = entry
        parts = ['helix', 'loop']
        value_names = ['bps', 'stacks', 'bphs']
        summary = Summary(pdb=pdb, chain=chain)

        for range_name in ['lr', 'sr']:
            for part1 in parts:
                for part2 in parts:
                    values = self.query(pdb, chain, part1, part2, range_name)
                    for name, value in zip(value_names, values):
                        full = '%s_%s_%s_%s' % (range_name, part1, part2, name)
                        setattr(summary, full, value)

        return summary
