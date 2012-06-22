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

from MLSqlAlchemyClasses import session, LoopSearch, LoopSearchQA, LoopPositions
from MotifAtlasBaseClass import MotifAtlasBaseClass


class LoopSearchesLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        MotifAtlasBaseClass._setup_matlab(self)

        self.loopSearchDir = self.config['locations']['loops_search_dir']
        self.precomputedData = self.config['locations']['loops_mat_files']
        self.loop_regex = '(IL|HL)_\w{4}_\d{3}'
        self.pdb_regex = '^[0-9A-Za-z]{4}$'
        self.update = True # determines whether to update existing values in the db

    def load_loop_annotations(self):
        """
        """
        # loop over directories
        for folder in os.listdir(self.precomputedData):
            if re.search(self.pdb_regex, folder):
                logging.info('Importing loop annotations from %s', folder)
            else:
                continue
            [outputFile, err_msg] = self.mlab.loadLoopAnnotations(os.path.join(self.precomputedData, folder), nout=2)
            if err_msg != '':
                MotifAtlasBaseClass._crash(self, err_msg)
            else:
                reader = csv.reader(open(outputFile), delimiter=',', quotechar='"')
                for row in reader:
                    (loop_id, position, nt_id, bulge, flanking) = row
                    existing = session.query(LoopPositions). \
                                       filter(LoopPositions.loop_id==loop_id). \
                                       filter(LoopPositions.position==position). \
                                       first()
                    if existing:
                        if self.update:
                            existing.flanking = int(flanking)
                            existing.bulge = int(bulge)
                            existing.nt_id = nt_id
                            session.merge(existing)
                    else:
                        session.add(LoopPositions(loop_id=loop_id,
                                                  position=position,
                                                  nt_id=nt_id,
                                                  flanking=int(flanking),
                                                  bulge=int(bulge)))
                session.commit()
                os.remove(outputFile) # delete temporary csv file

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

    def load_loop_search_qa_text_file(self, file):
        """
        """
        reader = csv.reader(open(file, 'r'))
        for row in reader:
            existing = session.query(LoopSearchQA). \
                               filter(LoopSearchQA.loop_id1==row[0]). \
                               filter(LoopSearchQA.loop_id2==row[1]). \
                               first()
            if existing:
                if self.update:
                    existing.status = int(row[2])
                    existing.message = row[3]
                    session.merge(existing)
            else:
                session.add(LoopSearchQA(loop_id1=row[0],
                                         loop_id2=row[1],
                                         status=int(row[2]),
                                         message=row[3]))
        session.commit()


def usage():
    print __doc__

def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)
    S = LoopSearchesLoader()

#     S.load_loop_searches()

#     S.load_loop_search_qa_text_file('/Users/anton/FR3D/MM_extraNTs.txt')

#     S.load_loop_search_qa_text_file('/Users/anton/FR3D/MM_symmetrize.txt')

    S.load_loop_annotations()

#     S.load_final_matching_matrix()


if __name__ == "__main__":
    main(sys.argv[1:])