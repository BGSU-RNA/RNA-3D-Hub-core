"""

About

"""

__author__ = 'Anton Petrov'

import pdb, sys, getopt, logging, datetime

from MLSqlAlchemyClasses import session, LoopQA, LoopRelease, LoopModifications
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

        [result, L, err_msg] = self.mlab.aLoopQualityAssurance(pdb_id, nout=3)

        if err_msg != '':
            MotifAtlasBaseClass._crash(self,err_msg)

        for i in xrange(L):
            session.add(LoopQA(id   = result[i].id,
                               code = result[i].status,
                               self_compl = result[i].self_compl,
                               release_id = release_id))
            if result[i].modres != '':
                session.add(LoopModifications(id = result[i].id,
                                              modification = result[i].modres,
                                              release_id   = release_id))
        session.commit()
        self.mark_pdb_as_analyzed(pdb_id,'qa')


def usage():
    print __doc__

def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)
    Q = LoopQualityChecker()

#     pdbs = ['1HLX','124D','2AW4']
    pdbs = ['1HLX']

    Q.check_loop_quality(pdbs)


if __name__ == "__main__":
    main(sys.argv[1:])