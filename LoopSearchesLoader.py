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
import csv

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
        self.update = True # determines whether to update existing values in the db


    def load_loop_searches(self):
        """
            directory structure: loopSearchDir filesep IL_1S72_001 filesep IL_1S72_001_IL_1J5E_001.mat
        """
        # loop over directories
        for folder in os.listdir(self.loopSearchDir):
            if re.search(self.loop_regex, folder):
                logging.info('Importing %s searches', folder)
            else:
                continue
            # run matlab to create a temporary csv file with results
            [outputFile, err_msg] = self.mlab.loadLoopSearchFile(os.path.join(self.loopSearchDir, folder), nout=2)
            if err_msg != '':
                MotifAtlasBaseClass._crash(self, err_msg)
            else:
                reader = csv.reader(open(outputFile), delimiter=',', quotechar='"')
                for row in reader:
                    (loop_id1, loop_id2, disc, nt_list1, nt_list2) = row
                    self._store_in_database(loop_id1, loop_id2, disc, nt_list1, nt_list2)
                os.remove(outputFile) # delete temporary csv file

            # read in No_candidates.txt if exists
            self._read_no_candidates_file(folder)

    def _store_in_database(self, loop_id1, loop_id2, disc=-1, nt_list1=None, nt_list2=None):
        """
        """
        existing = session.query(LoopSearch). \
                           filter(LoopSearch.loop_id1==loop_id1). \
                           filter(LoopSearch.loop_id2==loop_id2). \
                           first()
        if existing:
            if self.update:
                existing.disc = float(disc)
                existing.nt_list1 = nt_list1
                existing.nt_list2 = nt_list2
                session.merge(existing)
        else:
            session.add(LoopSearch(loop_id1=loop_id1,
                                   loop_id2=loop_id2,
                                   disc=float(disc),
                                   nt_list1=nt_list1,
                                   nt_list2=nt_list2))
        session.commit()

    def _read_no_candidates_file(self, folder):
        """
        """
        no_candidates_file = os.path.join(self.loopSearchDir, folder, 'No_candidates.txt')
        if os.path.exists(no_candidates_file):
            loop_id1 = folder
            loops = open(no_candidates_file, 'r').readlines()
            loops = [x.rstrip() for x in loops]
            for loop_id2 in loops:
                self._store_in_database(loop_id1, loop_id2)

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