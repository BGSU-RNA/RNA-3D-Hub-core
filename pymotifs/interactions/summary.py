from pymotifs import core

from pymotifs.models import ChainInfo
from pymotifs.models import PdbHelixLoopInteractionSummary as Summary


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


class Loader(core.SimpleLoader):

    def query(self, session, pdb):
            return session.query(Summary).filter_by(pdb=pdb)

    def chains(self, pdb):
        with self.session() as session:
            query = session.query(ChainInfo).filter_by(pdb_id=pdb)
            if not query.count():
                raise core.SkipPdb("No chains found for %s" % pdb)

            return [result.chain_name for result in query]

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

    def result(self, pdb, chain, element1, element2, range_type):
        query = self.build(element1, element2, range_type)
        with self.session() as session:
            results = session.execute(query, {'pdb': pdb, 'chain': chain})

            results = results.fetchone()
            if not results:
                raise core.SkipValue("Couldn't compute %s %s" %
                                     element1, element2)

            return results['bps'], results['stacks'], results['bphs']

    def data(self, pdb, **kwargs):
        parts = ['helix', 'loop']
        for chain in self.chains(pdb):
            summary = Summary(pdb=pdb, chain=chain)

            for range_name in ['lr', 'sr']:
                for part1 in parts:
                    for part2 in parts:
                        values = self.result(pdb, chain, part1, part2,
                                             range_name)

                        for name, value in values.items():
                            full = '_'.join([range_name, part1, part2, name])
                            setattr(summary, full, value)

        return summary
