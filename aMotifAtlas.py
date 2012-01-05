"""

About

"""

__author__ = 'Anton Petrov'

import sys, getopt, logging

from aMLSqlAlchemyClasses import *
from aDistancesAndCoordinatesLoader import DistancesAndCoordinatesLoader
from aPdbInfoLoader import PdbInfoLoader
from aLoopExtractor import LoopExtractor
from aLoopQualityChecker import LoopQualityChecker


def usage():
    print __doc__

def main(argv):
    """
    """
#     logging.basicConfig(level=logging.DEBUG)
    logfile = 'motifatlas.log'
    logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)
    logging.info('Initializing update')

#     pdbs = ['1EKA','1HLX']#,'1S72','2AVY']

    """get new pdb files, import descriptions into the database"""
    P = PdbInfoLoader()
    P.get_all_rna_pdbs()
#     P.update_rna_containing_pdbs()
    """import coordinates and distances into the database"""
    D = DistancesAndCoordinatesLoader()
    D.import_distances(P.pdbs[1:200], recalculate=False)
    D.import_coordinates(P.pdbs[1:200],recalculate=False)
    """extract all loops and import into the database"""
    E = LoopExtractor()
    E.extract_and_import_loops(P.pdbs[1:200])
    """do loop QA, import into the database"""
    Q = LoopQualityChecker()
    Q.check_loop_quality(P.pdbs[1:200])

    # compute new non-redundant lists, import into the database

    # get a list of loops to be clustered

    # cluster motifs, import into the database

    # annotate all pdb files with these clusters

    # on failure: stop, email

    # log from matlab, log from python

    # log with dates, clear filenames



    logging.info('SUCCESSFUL UPDATE')
    E.send_report('motifatlas.log')

if __name__ == "__main__":
    main(sys.argv[1:])