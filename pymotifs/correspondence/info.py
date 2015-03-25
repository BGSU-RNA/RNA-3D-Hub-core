from pymotifs import core
from pymotifs.models import CorrespondenceInfo as Info

from sqlalchemy.sql.expression import text


# This uses a cross join which is very expensive but should be ok here as there
# are relatively few unique rna sequences.
QUERY = """
select distinct
    if(E1.id < E2.id, E1.id, E2.id),
    if(E1.id < E2.id, E2.id, E1.id)
from exp_seq_info as E1
join exp_seq_info as E2
join exp_seq_chain_mapping as M1
on
    M1.exp_seq_id = E1.id
join exp_seq_chain_mapping as M2
on
    M2.exp_seq_id = E2.id
join chain_info as I1
on
    I1.id = M1.chain_id
join chain_info as I2
on
    I2.id = M2.chain_id
where
    I1.pdb_id in ({items})
    and I2.pdb_id in ({items})
    and (
    (
        E1.length < 36
        and E2.length = E1.length
    )
    or (
        E2.length <= 2 * E1.length
        and E2.length >= (E1.length / E2.length)
    )
    )
;
"""


class Loader(core.SimpleLoader):
    """A class to load up all pairs of experimental sequences that should have
    an alignment attempted. This does not work per structure as many other
    things do, but instead will compute all possible pairs and then work per
    pair, inserting or storing each pair as needed.
    """

    short_cutoff = 36
    allow_no_data = True
    mark = False

    def to_process(self, pdbs, **kwargs):
        """Here we transform the list of pdbs to a list of pairs of
        experimental sequences to try and align. This simplifies a lot of the
        logic later on, but does cause problems with marking stuff as
        processed.

        :pdbs: The list of pdbs to process.
        """

        items = ','.join(['"%s"' % pdb for pdb in pdbs])

        pairs = []
        with self.session() as session:
            query = session.execute(text(QUERY.format(items=items)))
            for result in query.fetchall():
                pairs.append((result[0], result[1]))

        return pairs

    def query(self, session, pair, **kwargs):
        return session.query(Info).\
            filter_by(exp_seq_id1=pair[0], exp_seq_id2=pair[1])

    def data(self, pair, **kwargs):
        return Info(exp_seq_id1=pair[0], exp_seq_id2=pair[1])
