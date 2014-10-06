"""

Program for importing loop quality assurance data into the RNA 3D Hub database.

Some loops should be disqualified because they have unresolved or missing
nucleotides. See matlab code for disqualification codes.

"""

__author__ = 'Anton Petrov'

import pdb
import sys
import getopt
import logging
import datetime
import csv
import os


from models import session, LoopQa, LoopReleases
from MotifAtlasBaseClass import MotifAtlasBaseClass
from PdbInfoLoader import PdbInfoLoader

logger = logging.getLogger(__name__)


class LoopQualityChecker(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)

    def check_loop_quality(self, pdbs):
        """
        """
        try:
            logger.info('Loop Quality Assurance')
            release = LoopReleases(mode=self.config['release_mode']['loops'])
            for pdb_id in pdbs:
                self.loop_qa(pdb_id, release.id)
            session.add(release)
            session.commit()
            logger.info('Loop QA complete')
            logger.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def loop_qa(self, pdb_id, release_id):
        """
        """
        logger.info('QA on %s', pdb_id)
        MotifAtlasBaseClass._setup_matlab(self)

        [ifn, err_msg] = self.mlab.aLoopQualityAssurance(pdb_id, nout=2)

        if err_msg != '':
            logger.warning('Error %s in pdb %s' % (err_msg, pdb_id))
        else:
            self.__import_qa_from_csv(ifn, release_id)
            self.mark_pdb_as_analyzed(pdb_id,'qa')

    def __import_qa_from_csv(self, ifn, release_id):
        """Reads the csv file, imports all distances, deletes the file when done
           to avoid stale data and free up disk space"""
        logger.info('Importing qa')
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        QA = []
        for i, row in enumerate(reader):
            modres = row[2]
            if modres == '': modres = None
            compl = row[4]
            if compl == '': compl = None

            QA.append(LoopQa(id     = row[0],
                             status = row[1],
                             modifications = modres,
                             nt_signature  = row[3],
                             complementary = compl,
                             release_id = release_id))
        os.remove(ifn)
        session.add_all(QA)
        session.commit()
        logger.info('%s loops checked and imported' % len(QA))


def usage():
    print __doc__

def main(argv):
    """
    """
    Q = LoopQualityChecker()
    Q.start_logging()

#     pdbs = ['1FG0']

    p = PdbInfoLoader()
    p.get_all_rna_pdbs()
    pdbs = p.pdbs

    Q.check_loop_quality(pdbs)


if __name__ == "__main__":
    main(sys.argv[1:])
