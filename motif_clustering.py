"""

About

"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging, datetime
from MLSqlAlchemyClasses import session, AllLoops
from MotifAtlasBaseClass import MotifAtlasBaseClass
# from NRSqlAlchemyClasses import


class ClusterMotifs(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.pdbs  = []
        self.loops = []

    def get_pdb_ids_for_clustering(self):
        """
        """
        self.pdbs.append('1J5E')

    def get_loops_for_clustering(self):
        """
        """
        self.loops = session.query(AllLoops).filter(AllLoops.pdb.in_(self.pdbs)).all()
        print self.loops


#     def extract_and_import_loops(self):
#         """Loops over `pdbs`, extracts and imports all loops"""
#         try:
#
#         except:
#             e = sys.exc_info()[1]
#             MotifAtlasBaseClass._crash(self,e)



def usage():
    print __doc__

def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)
    M = ClusterMotifs()

    M.get_pdb_ids_for_clustering()
    M.get_loops_for_clustering()


if __name__ == "__main__":
    main(sys.argv[1:])

