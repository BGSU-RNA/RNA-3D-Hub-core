"""
Read the nr_classes and nr_chains tables and write a new nr_class_rank table to replace nr_chains
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict


class class_rank_populate_loader(core.SimpleLoader):
    """
    Loader to store quality data for an input equivalence class
    in table nr_cqs.
    """

    dependencies = set([])

    merge_data = True

    """We allow for no data to be written when appropriate"""
    allow_no_data = True

    def to_process(self, pdbs, **kwargs):
        """
        Just return a list with one item
        """

        return [1]

    def query(self, session, nr_name):
        return session.query(mod.NrCqs.nr_name).\
            filter(mod.NrCqs.nr_name == nr_name)

    def data(self, nr_name, **kwargs):
        """
        Read the nr_classes and nr_chains tables.
        Keep just the most recent ranking.
        Write to the new nr_class_rank table.
        """

        with self.session() as session:
            # join the two tables, order by nr_class_id so that the most recent is last
            query = session.query(mod.NrClasses.name,mod.NrChains.ife_id,mod.NrChains.rank).\
                join(mod.NrChains, mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
                order_by(mod.NrClasses.nr_class_id)

            # map class_name and ife to rank
            name_ife_to_rank = {}
            for result in query:
                name = result.name
                ife = result.ife_id
                rank = result.rank
                name_ife = name + "&" + ife
                name_ife_to_rank[name_ife] = rank

            # sort by name, then rank, then ife
            name_rank_ife = []
            for name_ife, rank in name_ife_to_rank.items():
                name, ife = name_ife.split("&")
                first_ife = ife.split("+")[0]
                if len(first_ife.split("|")) == 3:
                    name_rank_ife.append((name, rank, ife))

            # write to database in the nicest order
            for name, rank, ife in sorted(name_rank_ife):
                print("%s %s %s" % (name, ife, rank))
                self.logger.info("%s %s %s" % (name, ife, rank))

                yield mod.NrClassRank(
                    nr_class_name = name,
                    ife_id = ife,
                    rank = rank)


