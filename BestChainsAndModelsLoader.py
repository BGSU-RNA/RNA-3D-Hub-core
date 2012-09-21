"""

Import information about redundancy within PDB files.

"""

__author__ = 'Anton Petrov'

import os
import csv
import sys
import logging

from MLSqlAlchemyClasses import session, PdbBestChainsAndModels, PdbAnalysisStatus
from MotifAtlasBaseClass import MotifAtlasBaseClass


class BestChainsAndModelsLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)

    def import_best_chains_and_models(self, pdbs, recalculate=False):
        """
        """
        try:
            logging.info('Importing best chains and models')
            if not recalculate:
                recalculate = self.config['recalculate']['best_chains_and_models']
            if recalculate:
                pdb_list = pdbs
                self.__delete_old_data(pdbs)
            else:
                pdb_list = self.filter_out_analyzed_pdbs(pdbs,'best_chains_and_models')

            if pdb_list:
                MotifAtlasBaseClass._setup_matlab(self)

            for pdb_file in pdb_list:
                logging.info('Running matlab on %s', pdb_file)
                # 'ABC', '1,2', ''
                best_chains, best_models, err_msg = self.mlab.loadBestChainsAndModels(pdb_file, nout=3)
                best_chains = ','.join(list(best_chains))

                if err_msg == '':
                    self.__import_into_db(pdb_file, best_chains, best_models)
                else:
                    MotifAtlasBaseClass._crash(self,err_msg)

                self.mark_pdb_as_analyzed(pdb_file,'best_chains_and_models')

            logging.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def __import_into_db(self, pdb_id, best_chains, best_models):
        """
        """
        try:
            session.add(PdbBestChainsAndModels(pdb_id=pdb_id,
                                               best_chains=best_chains,
                                               best_models=best_models))
            logging.info('Import complete for %s' % pdb_id)
        except:
            # do nothing in case of duplicate entry
            pass
        session.commit()

    def __delete_old_data(self, pdbs):
        """When recalculate=True, delete all data from pdb_best_chains_and_models
        and set pdb_analysis_status to None in the best_chains_and_models column"""
        session.query(PdbBestChainsAndModels). \
                filter(PdbBestChainsAndModels.pdb_id.in_(pdbs)). \
                delete(synchronize_session='fetch')
        for statusObj in session.query(PdbAnalysisStatus). \
                                 filter(PdbAnalysisStatus.id.in_(pdbs)). \
                                 all():
            statusObj.best_chains_and_models = None
            session.merge(statusObj)
        session.commit()


def usage():
    print __doc__

def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    B = BestChainsAndModelsLoader()

    from PdbInfoLoader import PdbInfoLoader
    P = PdbInfoLoader()
    P.get_all_rna_pdbs()

    B.import_best_chains_and_models(P.pdbs)
#     B.import_best_chains_and_models(P.pdbs, recalculate=True)


if __name__ == "__main__":
    main(sys.argv[1:])
