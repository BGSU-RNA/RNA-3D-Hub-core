from sqlalchemy.orm import aliased
from pymotifs import core
from pymotifs import models as mod

from pymotifs.chain_chain.comparison import Loader as ChainComparisons

class Exporter(core.Exporter):
    '''Export discrepancy values of every possible ife id combination.
    '''

    # General setup

    dependencies = set([ChainComparisons])

    def filename(self, *args, **kwargs):
        """This will always return the configured path at
        locations.interactions_gz and is where the interaction export will be
        written.

        Returns
        -------
        filename : str
            The path to write to.
        """
        return None

    def to_process(self, pdbs, **kwargs):

        if len(pdbs) < 500:
            raise core.Skip("Too few pdb files being processed to need to write discrepancies")

        # return a list of length 1, so it only writes out once
        return ["pretend_PDB_ID"]


    def process(self, pdb, **kwargs):
        """ override the process method to simply load data,
        which takes care of the writing, and nothing else
        needs to be returned or written
        """

        data = self.data(pdb)

    def data(self, pdb, **kwargs):

        '''This method dumps the discrepancy for every possible combination of ife ids into a file for later import.
        The data is pickled for easy reading and writing in python.

            :results: list of lists of the the two ife ids and their corresponding discrepancy.
        '''

        with self.session() as session:
            self.logger.info("ife_discrepancy_dump: Inside data retrieval routine")
            self.logger.info("ife_discrepancy_dump: Setting aliases")

            IC1 = aliased(mod.IfeChains)
            IC2 = aliased(mod.IfeChains)
            CCS = mod.ChainChainSimilarity

            self.logger.info("ife_discrepancy_dump: Building Query")

            query = session.query(CCS.discrepancy,
                                  IC1.ife_id.label('ife1'),
                                  IC2.ife_id.label('ife2')).\
                                  join(IC1, IC1.chain_id == CCS.chain_id_1).\
                                  join(IC2, IC2.chain_id == CCS.chain_id_2).\
                                  filter(IC1.index == 0).\
                                  filter(IC2.index == 0).\
                                  filter(IC1.chain_id != IC2.chain_id).\
                                  order_by(IC1.ife_id, IC2.ife_id).\
                                  distinct()
#            query.all()

            self.logger.info("ife_discrepancy_dump: Query Built")

            count = query.count()
            if not count:
                self.logger.warning("No discrepancies found.")
            else:
                self.logger.info("Found %s discrepancies.", count)

            # create a list of lists made of the query results

            results = [[r.ife1, r.ife2, r.discrepancy] for r in query]

            # write into a file
            with open('/var/www/html/discrepancy/IFEdiscrepancy.txt', 'w') as f:
                self.logger.info("ife_discrepancy_dump: File Open")
                for r in results:
                    f.write('%s\t%s\t%s\n' % (r[0], r[1], r[2]))

        return []
