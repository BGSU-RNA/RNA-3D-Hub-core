"""

About

"""

__author__ = 'Anton Petrov'

import pdb, sys, getopt, logging, datetime, csv, os

from models import session, LoopQA, LoopRelease
from MotifAtlasBaseClass import MotifAtlasBaseClass


class LoopQualityChecker(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)

    def check_loop_quality(self, pdbs):
        """
        """
        try:
            logging.info('Loop Quality Assurance')
            release = LoopRelease(mode=self.config['release_mode']['loops'])
            for pdb_id in pdbs:
                self.loop_qa(pdb_id, release.id)
            session.add(release)
            session.commit()
            logging.info('Loop QA complete')
            logging.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def loop_qa(self, pdb_id, release_id):
        """
        """
        logging.info('QA on %s', pdb_id)
        MotifAtlasBaseClass._setup_matlab(self)

        [ifn, err_msg] = self.mlab.aLoopQualityAssurance(pdb_id, nout=2)

        if err_msg != '':
            MotifAtlasBaseClass._crash(self,err_msg)
        else:
            self.__import_qa_from_csv(ifn, release_id)
            self.mark_pdb_as_analyzed(pdb_id,'qa')

    def __import_qa_from_csv(self, ifn, release_id):
        """Reads the csv file, imports all distances, deletes the file when done
           to avoid stale data and free up disk space"""
        logging.info('Importing qa')
        reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
        QA = []
        for i, row in enumerate(reader):
            modres = row[2]
            if modres == '': modres = None
            compl = row[4]
            if compl == '': compl = None

            QA.append(LoopQA(id     = row[0],
                             status = row[1],
                             modifications = modres,
                             nt_signature  = row[3],
                             complementary = compl,
                             release_id = release_id))
        os.remove(ifn)
        session.add_all(QA)
        session.commit()
        logging.info('Csv file successfully imported')


def usage():
    print __doc__

def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)
    Q = LoopQualityChecker()

#     pdbs = ['1HLX','124D','2AW4']
#     pdbs = ['1HLX']
    pdbs = ['1ASY']

    Q.check_loop_quality(pdbs)


if __name__ == "__main__":
    main(sys.argv[1:])