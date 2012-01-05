"""

About

"""

__author__ = 'Anton Petrov'

import os, csv, pdb, sys, getopt, logging
from email.mime.text import MIMEText
from datetime import date

from aMLSqlAlchemyClasses import *
from aDistancesAndCoordinatesLoader import DistancesAndCoordinatesLoader
# from aLoopLoader    import LoopLoader
from aPdbInfoLoader import PdbInfoLoader
from aLoopExtractor import LoopExtractor
from aLoopQualityChecker import LoopQualityChecker

def usage():
    print __doc__

def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)
#     logging.basicConfig(filename='motifatlas.log', filemode='w', level=logging.DEBUG)
    logging.info('Initializing update')

#     pdbs = ['1EKA','1HLX']#,'1S72','2AVY']


    """get new pdb files, import descriptions into the database"""
    P = PdbInfoLoader()
    P.get_all_rna_pdbs()
#     P.update_rna_containing_pdbs()
    """import coordinates and distances into the database"""
    D = DistancesAndCoordinatesLoader()
    D.import_distances(P.pdbs[1:5], recalculate=False)
    D.import_coordinates(P.pdbs[1:5],recalculate=False)
    """extract all loops and import into the database"""
    E = LoopExtractor()
    E.extract_and_import_loops(P.pdbs[1:5])
    """do loop QA, import into the database"""
    Q = LoopQualityChecker()
    Q.check_loop_quality(P.pdbs[1:5])

# compute new non-redundant lists, import into the database



# get a list of loops to be clustered

# cluster motifs, import into the database

# annotate all pdb files with these clusters

# on failure: stop, email

# log from matlab, log from python

# log with dates, clear filenames



    logging.info('SUCCESSFUL UPDATE')
#     E.send_report()

if __name__ == "__main__":
    main(sys.argv[1:])


#     try:
#         opts, args = getopt.getopt(argv, "", ['help'])
#     except getopt.GetoptError:
#         usage()
#         sys.exit(2)
#     for opt, arg in opts:
#         pass
#         if   opt == '-d':
#             M.import_distances()
# #         elif opt == '-q':
# #             L.import_loop_qa()
# #         elif opt == '-c':
# #             L.import_coordinates()
# #         elif opt == '-l':
# #             L.import_all_loops()
# #         elif opt == '-m':
# #             L.import_loop_modifications()
# #         elif opt == '-r':
# #             L.remove_loop_release(arg)
#         elif opt in ('-h', '--help'):
#             usage()
#             sys.exit()