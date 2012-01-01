"""

About

"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging
from aMLSqlAlchemyClasses import session, AllLoops


class LoopExtractor:
    """
    """
    def __init__(self):
        self.mlab = False
        self.save_location = '/Users/anton/FR3D/MotifAtlas/PrecomputedData'
        if not os.path.exists(self.save_location):
            os.makedirs(self.save_location)

    def _setup_matlab(self):
        if self.mlab:
            logging.info('Matlab already running')
            return
        logging.info('Starting up matlab')
        from mlabwrap import mlab
        self.mlab = mlab
        self.mlab.setup()
        logging.info('Matlab started')

    def __crash(self, msg=None):
        if msg:
            logging.warning(msg)
        logging.critical('PROGRAM %s CRASHED', __name__)
        session.rollback()
        sys.exit(2)

    def get_loops(self, pdb_id, loop_type):
        """
        """
        try:
            if self._check_if_done(pdb_id, loop_type):
                logging.info('%s already extracted from pdb %s', loop_type, pdb_id)
                return

            self._setup_matlab()
            [Loops, l, err_msg] = self.mlab.aGetLoops(pdb_id,loop_type,nout=3)

            if err_msg != '':
                self.__crash(err_msg)

            for i in xrange(l):
                loop_id = self._get_loop_id(Loops[i].AllLoops_table.full_id,
                                            pdb_id, loop_type)
                Loops[i].Filename = loop_id
                session.add(AllLoops(id            = loop_id,
                                     type          = loop_type,
                                     pdb           = pdb_id,
                                     sequential_id = loop_id[-3:],
                                     length        = int(Loops[i].NumNT[0][0]),
                                     seq           = Loops[i].AllLoops_table.seq,
                                     r_seq         = Loops[i].AllLoops_table.r_seq,
                                     nwc_seq       = Loops[i].AllLoops_table.nwc,
                                     r_nwc_seq     = Loops[i].AllLoops_table.r_nwc,
                                     pdb_file      = Loops[i].PDBFilename,
                                     nt_ids        = Loops[i].AllLoops_table.full_id,
                                     loop_name     = Loops[i].AllLoops_table.loop_name))

            session.commit()
            logging.info('%s from %s successfully imported', loop_type, pdb_id)
#             self.mlab.aSaveLoops(Loops, ids, self.save_location)
        except:
            e = sys.exc_info()[1]
            logging.critical('matlab_import_distances CRASHED')
            self.__crash(e)

    def _check_if_done(self, pdb_id, loop_type):
        """
        """
        if session.query(AllLoops). \
                   filter(AllLoops.pdb==pdb_id). \
                   filter(AllLoops.type==loop_type).first():
            return True
        else:
            return False

    def _get_loop_id(self, full_id, pdb_id, loop_type):
        """
        """
        L = session.query(AllLoops).filter(AllLoops.nt_ids==full_id).first()
        if L:
            logging.info('Full_id %s matched %s', full_id, L.id)
            return L.id
        else:
            "count the loops already in the db"""
            seq_id = session.query(AllLoops). \
                             filter(AllLoops.pdb==pdb_id). \
                             filter(AllLoops.type==loop_type).count()
            """format example: IL_1S72_001"""
            id = '_'.join([loop_type, pdb_id, str(seq_id+1).rjust(3,'0')])
            logging.info('Created new loop id %s', id)
            return id


def usage():
    print __doc__

def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)
    E = LoopExtractor()

#     pdbs = ['1S72','124D','2AW4']
    pdbs = ['1HLX']#,'2AW4']

    for pdb_id in pdbs:
        E.get_loops(pdb_id,'IL')
        E.get_loops(pdb_id,'HL')
#         E.get_loops(pdb_id,'J3')


if __name__ == "__main__":
    main(sys.argv[1:])