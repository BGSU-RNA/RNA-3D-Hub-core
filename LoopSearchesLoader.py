"""

About

"""

__author__ = 'Anton Petrov'

import pdb
import sys
import getopt
import logging
import os
import re

from MLSqlAlchemyClasses import session, LoopSearch
from MotifAtlasBaseClass import MotifAtlasBaseClass


class LoopSearchesLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        MotifAtlasBaseClass._setup_matlab(self)

        self.loopSearchDir = self.config['locations']['loops_search_dir']
        self.loop_regex = '(IL|HL)_\w{4}_\d{3}'
        self.matfile_regex = '(%s)_(%s)\.mat' % (self.loop_regex, self.loop_regex)
        self.loop_id1 = self.loop_id2 = ''
        self.update = True # determines whether to update existing values in the db


    def load_loop_searches(self):
        """
            directory structure: loopSearchDir filesep IL_1S72_001 filesep IL_1S72_001_IL_1J5E_001.mat
        """
        # loop over directories
        for dir in os.listdir(self.loopSearchDir):
            if re.search(self.loop_regex, dir):
                self.loop_id1 = dir
                logging.info('Importing %s searches', self.loop_id1)
            else:
                continue
            # loop over files in each directory
            for file in os.listdir( os.path.join(self.loopSearchDir, self.loop_id1) ):
                fullFile = os.path.join(self.loopSearchDir, dir, file)
                match = re.search(self.matfile_regex, file)
                if match:
                    self.loop_id2 = match.group(3)
                    [disc, nt_list1, nt_list2, err_msg] = self.mlab.loadLoopSearchFile(fullFile, nout=4)
                    if err_msg != '':
                        MotifAtlasBaseClass._crash(self, err_msg)
                    else:
                        self._store_in_database(disc[0][0], nt_list1, nt_list2)
                elif file == 'No_candidates.txt':
                    self._read_no_candidates_file(fullFile)
                else:
                    continue

    def _store_in_database(self, disc=-1, nt_list1=None, nt_list2=None):
        """
        """
        existing = session.query(LoopSearch). \
                           filter(LoopSearch.loop_id1==self.loop_id1). \
                           filter(LoopSearch.loop_id2==self.loop_id2). \
                           first()
        if existing:
            if self.update:
                existing.disc = disc
                existing.nt_list1 = nt_list1
                existing.nt_list2 = nt_list2
                session.merge(existing)
        else:
            session.add(LoopSearch(loop_id1=self.loop_id1,
                                   loop_id2=self.loop_id2,
                                   disc=disc,
                                   nt_list1=nt_list1,
                                   nt_list2=nt_list2))
        session.commit()

    def _read_no_candidates_file(self, fullFile):
        """
        """
        loops = open(fullFile , 'r').readlines()
        loops = [x.rstrip() for x in loops]
        for self.loop_id2 in loops:
            self._store_in_database()

    def load_final_matching_matrix(self):
        """
        """
        pass


def usage():
    print __doc__

def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)
    S = LoopSearchesLoader()

    S.load_loop_searches()
#     S.load_final_matching_matrix()


if __name__ == "__main__":
    main(sys.argv[1:])