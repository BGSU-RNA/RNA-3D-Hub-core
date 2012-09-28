"""

Import information about redundancy between nucleotides within a PDB file


"""

__author__ = 'Anton Petrov'

import os
import csv
import sys
import logging

from models import session, RedundantNucleotide, PdbAnalysisStatus
from MotifAtlasBaseClass import MotifAtlasBaseClass


class RedundantNucleotidesLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)

    def import_redundant_nucleotides(self, pdbs, recalculate=False):
        """
        """
        try:
            logging.info('Importing redundant nucleotides')
            if not recalculate:
                recalculate = self.config['recalculate']['redundant_nts']
            if recalculate:
                pdb_list = pdbs
                self.__delete_old_data(pdbs)
            else:
                pdb_list = self.filter_out_analyzed_pdbs(pdbs,'redundant_nts')

            if pdb_list:
                MotifAtlasBaseClass._setup_matlab(self)

            for pdb_file in pdb_list:
                logging.info('Running matlab on %s', pdb_file)
                ifn, err_msg = self.mlab.loadRedundantNucleotides(pdb_file, nout=2)
                if err_msg == '':
                    self.__import_temporary_file(ifn, pdb_file)
                else:
                    MotifAtlasBaseClass._crash(self,err_msg)

                self.mark_pdb_as_analyzed(pdb_file,'redundant_nts')

            logging.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def __import_temporary_file(self, ifn, pdb_file):
        """
        """
        logging.info('Importing redundant nucleotides')
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        for i, row in enumerate(reader):
            R = RedundantNucleotide(nt_id1=row[0],
                                        nt_id2=row[1],
                                        pdb_id=pdb_file)
            try:
                session.add(R)
            except:
                pass
        session.commit()
        os.remove(ifn)
        logging.info('File successfully imported')

    def __delete_old_data(self, pdbs):
        """When recalculate=True, delete all data from pdb_redundant_nucleotides
        and set pdb_analysis_status to None in the redundant_nts column"""
        session.query(RedundantNucleotide). \
                filter(RedundantNucleotide.pdb_id.in_(pdbs)). \
                delete(synchronize_session='fetch')
        for statusObj in session.query(PdbAnalysisStatus). \
                                 filter(PdbAnalysisStatus.id.in_(pdbs)). \
                                 all():
            statusObj.redundant_nts = None
            session.merge(statusObj)
        session.commit()


def usage():
    print __doc__

def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    R = RedundantNucleotidesLoader()

#     pdbs = ['1J5E', '1KOG']

    from PdbInfoLoader import PdbInfoLoader
    P = PdbInfoLoader()
    P.get_all_rna_pdbs()

    R.import_redundant_nucleotides(P.pdbs)
#     R.import_redundant_nucleotides(pdbs, recalculate=True)


if __name__ == "__main__":
    main(sys.argv[1:])
