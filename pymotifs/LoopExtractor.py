"""

Python module for extracting loops from RNA 3D structures.

"""

__author__ = 'Anton Petrov'

import os
import csv
import pdb
import sys
import getopt
import logging
import datetime


from models import session, AllLoops
from MotifAtlasBaseClass import MotifAtlasBaseClass

logger = logging.getLogger(__name__)


class LoopExtractor(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.loop_types = ['IL','HL','J3']

    def extract_and_import_loops(self, pdbs, recalculate=None):
        """Loops over `pdbs`, extracts and imports all loops"""
        try:
            for loop_type in self.loop_types:
                logger.info('Extracting %s' % loop_type)
                if recalculate is None:
                    recalculate = self.config['recalculate'][loop_type]
                if recalculate:
                    pdb_list = pdbs[:]
                else:
                    pdb_list = self.filter_out_analyzed_pdbs(pdbs, loop_type)
                for pdb_id in pdb_list:
                    logger.info('Extracting %s from %s', loop_type, pdb_id)
                    (Loops,l) = self.extract_loops(pdb_id, loop_type)
                    self.import_loops(Loops, l, pdb_id, loop_type)
                    logger.info('%s', '='*40)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def extract_loops(self, pdb_id, loop_type):
        """
        """
        try:
            MotifAtlasBaseClass._setup_matlab(self)
            """Loops - array of FR3D File structures. l - its length"""
            [Loops, l, err_msg] = self.mlab.extractLoops(pdb_id, loop_type, nout=3)

            if err_msg != '':
                MotifAtlasBaseClass._crash(self,err_msg)

            if Loops == 0:
                logger.info('No %s in %s', loop_type, pdb_id)
                return (0, 0)
            else:
                logger.info('Found %i loops', l)
                return (Loops, l)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def import_loops(self, Loops, l, pdb_id, loop_type):
        """
        """
        try:
            if Loops == 0:
                self.mark_pdb_as_analyzed(pdb_id, loop_type)
                return
            for i in xrange(l):
                loop_id = self._get_loop_id(Loops[i].AllLoops_table.full_id,
                                            pdb_id, loop_type)
                Loops[i].Filename = loop_id
                session.merge(
                    AllLoops(id            = loop_id,
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
            self.save_mat_files(Loops)
            self.mark_pdb_as_analyzed(pdb_id, loop_type)
            logger.info('%s from %s successfully imported', loop_type, pdb_id)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def save_mat_files(self,Loops):
        """Pass the Loops structure array back to matlab so that it can
        save the .mat files in the specified location."""
        MotifAtlasBaseClass._setup_matlab(self)
        [status, err_msg] = self.mlab.aSaveLoops(Loops,
                                                 self.config['locations']['loops_mat_files'],
                                                 nout=2)
        if status == 0:
            logger.info('mat files saved')
        else:
            MotifAtlasBaseClass._crash(self,err_msg)

    def _get_loop_id(self, full_id, pdb_id, loop_type):
        """returns a loop id"""
        L = session.query(AllLoops).filter(AllLoops.nt_ids==full_id).first()
        if L:
            logger.info('Full_id %s matched %s', full_id, L.id)
            return L.id
        else:
            "count the loops already in the db"""
            seq_id = session.query(AllLoops). \
                             filter(AllLoops.pdb==pdb_id). \
                             filter(AllLoops.type==loop_type).count()
            """format example: IL_1S72_001"""
            id = '_'.join([loop_type, pdb_id, str(seq_id+1).rjust(3,'0')])
            logger.info('Created new loop id %s', id)
            return id



def usage():
    print __doc__

def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)
    E = LoopExtractor()

#     pdbs = ['1HLX','124D','2AW4']
    pdbs = ['1HLX','2AW4']

    E.extract_and_import_loops(pdbs)# , recalculate=False)


if __name__ == "__main__":
    main(sys.argv[1:])


# SELECT * FROM loops_all AS t1
# LEFT JOIN motifversions.`loops_all` AS t2
# ON t1.id=t2.id
# WHERE t1.type!=t2.type
# OR t2.pdb!=t1.pdb
# OR t1.`sequential_id`!=t2.`sequential_id`
# OR t1.length!=t2.length
# OR t1.seq!=t2.seq
# OR t1.`r_nwc_seq`!=t2.`r_nwc_seq`
# OR t1.`r_seq`!=t2.`r_seq`
# OR t1.`nt_ids`!=t2.`nt_ids`
# OR t1.`loop_name` !=t2.`loop_name`
# OR t1.`pdb_file`!=t2.`pdb_file`
# OR t1.`nwc_seq`!=t2.`nwc_seq` #junctions are different;
