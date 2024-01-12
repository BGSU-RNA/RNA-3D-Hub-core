"""
Read the nr_classes and nr_chains tables and write a new nr_class_rank table to replace nr_chains
"""

from pymotifs import core
from pymotifs import models as mod
#from pymotifs.utils import row2dict


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
        Stratify by resolution
        """

        return ['all']
        return ['20.0']

        # all has too many to yield all at once
        return ['all_0','all_1','all_2','all_3','all_4','all_5','all_6','all_7','all_8','all_9']

        # return ['1.5','2.0','2.5','3.0','3.5','4.0','20.0','all']

    def query(self, session, nr_name):
        return session.query(mod.NrCqs.nr_name).\
            filter(mod.NrCqs.nr_name == nr_name)

    def data(self, resolution, **kwargs):
        """
        Read the nr_classes and nr_chains tables.
        Keep the most recent ranking, but keep the earliest nr_class_id
        Write to the new nr_class_rank table.
        """

        with self.session() as session:
            # join the two tables, order by nr_class_id so that the most recent is last
            query = session.query(mod.NrClasses.name,mod.NrClasses.nr_class_id,mod.NrChains.ife_id,mod.NrChains.rank).\
                join(mod.NrChains, mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
                filter(mod.NrClasses.name.like('NR_' + resolution + '%%')).\
                order_by(mod.NrClasses.nr_class_id)

            print('Query is done, processing rows')

            # map class_name and ife to rank
            name_ife_to_rank = {}
            name_ife_to_nr_class_id = {}
            c = 0
            for result in query:
                c += 1
                if c % 100000 == 0:
                    print('Processed %d rows from resolution %s' % (c,resolution))
                name = result.name
                ife = result.ife_id
                rank = result.rank
                name_ife = name + "&" + ife

                # overwrite previous rank, to use the most recently computed rank
                name_ife_to_rank[name_ife] = rank

                if not name_ife in name_ife_to_nr_class_id:
                    # only write once, presumably the first time it occurred
                    name_ife_to_nr_class_id[name_ife] = result.nr_class_id

        print('Done with query and saving')

        # sort by name, then rank, then ife
        name_rank_ife = []
        for name_ife, rank in name_ife_to_rank.items():
            id = name_ife_to_nr_class_id[name_ife]
            name, ife = name_ife.split("&")
            first_ife = ife.split("+")[0]

            # already done on rnatest, so don't do this again
            if len(first_ife.split("|")) == 3:
                name_rank_ife.append((name, id, rank, ife))

            # originally we did not do this on rnatest
            # release before 1.xx have no model numbers
            if len(first_ife.split("|")) == 2:
                name_rank_ife.append((name, id, rank, ife))

        # write to database in the nicest order
        c = 0
        for name, id, rank, ife in sorted(name_rank_ife):
            c += 1
            if c % 1000 == 0:
                print("%s %s %s %s" % (name, id, ife, rank))
            self.logger.info("%s %s %s %s" % (name, id, ife, rank))

            yield mod.NrClassRank(
                nr_class_name = name,
                nr_class_id = id,
                ife_id = ife,
                rank = rank)


