from pymotifs import core

from pymotifs.models import UnitInfo
from pymotifs.models import UnitPairsInteractions
from pymotifs.models import PdbHelixLoopInteractionSummary as Summary
from pymotifs.interactions.pairwise import Loader as InterLoader


QUERY_TEMPLATE = '''
SELECT
    count(if(f_lwbp regexp '^[ct]' {unique}, 1, NULL)) as 'bps',
    count(if(f_stacks regexp '^s' {unique}, 1, NULL)) as 'stacks',
    count(if(f_bphs regexp '^[1-9]'
            or (f_bphs = '0BPh' and unit1_id != unit2_id), 1, NULL)) as 'bphs'
FROM unit_pairs_interactions
JOIN {table1} as N1
ON
    unit1_id = N1.id
JOIN {table2} as N2
ON
    unit2_id = N2.id
WHERE
    pdb_id = :pdb
    AND N1.chain = N2.chain
    AND N1.chain = :chain
    AND (
        (f_lwbp regexp '^[ct]' {unique}) or
        (f_stacks regexp '^s' {unique}) or
        (f_bphs regexp '^[1-9]') or
        (f_bphs = '0BPh' and unit1_id != unit2_id)
    )
    AND f_crossing {operator} 4
;
'''


class Loader(core.SimpleLoader):
    dependencies = set([InterLoader])

    def has_data(self, pdb, **kwargs):                  # Check if we have already stored data for this pdb file in the database
        return True                                     # not sure why we have True values here, 

    def query(self, session, pdb):                              # this function does similar work with the has_data function
        return session.query(Summary).filter_by(pdb_id=pdb)

    def table(self, element):                               # key words changes
        if element == 'helix':
            return 'nt_helix_view'
        if element == 'loop':
            return 'nt_loop_view'
        raise ValueError("Unknown element type %s" % element)  

    def build(self, element1, element2, range_type):    # 
        if range_type == 'sr':
            operator = '<'
        elif range_type == 'lr':
            operator = '>='
        else:
            raise ValueError("Unknown range type: %s" % range_type)

        unique = ''
        if element1 == element2:
            unique = 'and unit1_id < unit2_id'

        return QUERY_TEMPLATE.format(unique=unique, operator=operator,
                                     table1=self.table(element1),
                                     table2=self.table(element2))

    def summary_query(self, pdb, chain, element1, element2, range_type):
        query = self.build(element1, element2, range_type)
        with self.session() as session:
            results = session.execute(query, {'pdb': pdb, 'chain': chain})

            results = results.fetchone()
            if not results:
                raise core.Skip("Couldn't compute %s %s", element1, element2)

            return results['bps'], results['stacks'], results['bphs']

    def chains(self, pdb):
        """Get all chains that have interactions. This only gets chains where
        the first interaction is in the chain.

        :pdb: The pdb id to get interactions for.
        :returns: A list of all chains for the given pdb file.
        """

        with self.session() as session:
            query = session.query(UnitInfo.chain).\
                join(UnitPairsInteractions,
                     UnitPairsInteractions.unit_id_1 == UnitInfo.unit_id).\
                filter(UnitInfo.pdb_id == pdb).\
                distinct()
            return [result.chain for result in query]

    def data(self, pdb, **kwargs):
        parts = ['helix', 'loop']
        value_names = ['bps', 'stacks', 'bphs']
        for chain in self.chains(pdb):
            summary = Summary(pdb_id=pdb, chain=chain)                      # I did not find where the Summary function was defined

            for range_name in ['lr', 'sr']:
                for part1 in parts:
                    for part2 in parts:
                        values = self.summary_query(pdb, chain, part1, part2,
                                                    range_name)
                        for name, value in zip(value_names, values):
                            full = '%s_%s_%s_%s' % (range_name, part1, part2,
                                                    name)
                            setattr(summary, full, value)
            yield summary
