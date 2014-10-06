"""

Program for loading pairwise interactions produced by FR3D into the RNA 3D Hub
database.

"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging
from models import session
from models import PdbPairwiseInteractions as PairwiseInteractions
from MotifAtlasBaseClass import MotifAtlasBaseClass

logger = logging.getLogger(__name__)


class PairwiseInteractionsLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.success = False

    def import_interactions(self, pdbs, recalculate=False):
        """Determines what files need to be analyzed, deletes stored data if
           necessary, loops over the pdbs, runs matlab on each of them
           independently, matlab generates a temporary csv file, it's imported
           and immediately deleted."""
        try:
            logger.info('Inside import_interactions')
            if not recalculate:
                recalculate = self.config['recalculate']['interactions']
            if recalculate:
                pdb_list = pdbs
                self.__delete_interactions(pdbs)
            else:
                pdb_list = self.filter_out_analyzed_pdbs(pdbs,'interactions')

            if pdb_list:
                MotifAtlasBaseClass._setup_matlab(self)

            for pdb_file in pdb_list:
                logger.info('Running matlab on %s', pdb_file)
                ifn, status, err_msg = self.mlab.loadInteractions(pdb_file,nout=3)
                status = status[0][0]
                if status == 0:
                    self.__import_interactions_from_csv(ifn, pdb_file)
                elif status == 2: # no nucleotides in the pdb file
                    logger.info('Pdb file %s has no nucleotides', pdb_file)
                else:
                    logger.warning('Matlab error code %i when analyzing %s',
                                     status, pdb_file)
                    MotifAtlasBaseClass._crash(self,err_msg)

                self.mark_pdb_as_analyzed(pdb_file,'interactions')
            self.success = True
            logger.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def __import_interactions_from_csv(self, ifn, pdb_file):
        """Reads the csv file, imports all interactions, deletes the file when
           done to avoid stale data and free up disk space"""
        logger.info('Importing interactions')
        commit_every = 1000
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for i,row in enumerate(reader):
            I = PairwiseInteractions(iPdbSig=row[0], jPdbSig=row[1])
            I.pdb_id = pdb_file
            I.f_crossing = int(row[3])
            interaction = row[2].strip()
            if interaction[0] == 's' or interaction[0:2] == 'ns':
                I.f_stacks = interaction
            elif interaction[-2:] == 'BR':
                I.f_brbs = interaction
            elif interaction[-3:] == 'BPh':
                I.f_bphs = interaction
            else:
                I.f_lwbp = interaction

            session.merge(I)
            """Since the files can be huge, it's unfeasible to store all
            objects in memory, have to commit regularly"""
            if i % commit_every == 0:
                session.commit()
        session.commit()
        os.remove(ifn)
        logger.info('Csv file successfully imported')

    def __delete_interactions(self, pdb_list):
        """recalculate=True, so delete what's already in the database"""
        logger.info('Deleting existing records %s', ','.join(pdb_list))
        for pdb_file in pdb_list:
            session.query(PairwiseInteractions). \
                    filter(PairwiseInteractions.pdb_id.in_(pdb_file+'%')). \
                    delete(synchronize_session=False)

def usage():
    print __doc__

def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    pdbs = ['3V11', '1J5E']

    from PdbInfoLoader import PdbInfoLoader
    p = PdbInfoLoader()
    p.get_all_rna_pdbs()
    pdbs = p.pdbs

    A = PairwiseInteractionsLoader()

    for pdb_id in pdbs:
        A.import_interactions([pdb_id], recalculate=True)



if __name__ == "__main__":
    main(sys.argv[1:])
