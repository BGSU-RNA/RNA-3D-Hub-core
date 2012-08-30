"""

Main entry point for launching a Motif Atlas update.

"""

__author__ = 'Anton Petrov'

import sys
import logging

from DistancesAndCoordinatesLoader import DistancesAndCoordinatesLoader
from PdbInfoLoader import PdbInfoLoader
from LoopExtractor import LoopExtractor
from LoopQualityChecker import LoopQualityChecker
from PairwiseInteractionsLoader import PairwiseInteractionsLoader
from RedundantNucleotidesLoader import RedundantNucleotidesLoader
from BestChainsAndModelsLoader import BestChainsAndModelsLoader
from MotifAtlasBaseClass import MotifAtlasBaseClass


def main(argv):
    """set up logging"""
    m = MotifAtlasBaseClass()
    m.start_logging()
    logging.info('Initializing update')

#     pdbs = ['1EKA','1HLX']#,'1S72','2AVY']

    """get new pdb files, import descriptions into the database"""
    p = PdbInfoLoader()
    p.get_all_rna_pdbs()
#     p.update_rna_containing_pdbs()
    p.check_obsolete_structures()
    """extract all loops and import into the database"""
    e = LoopExtractor()
    e.extract_and_import_loops(p.pdbs)
    """do loop QA, import into the database.
    This will always create a new loop release."""
#     q = LoopQualityChecker()
#     q.check_loop_quality(p.pdbs)
    """import pairwise interactions annotated by FR3D"""
    i = PairwiseInteractionsLoader()
    i.import_interactions(p.pdbs)
    """import coordinates and distances into the database"""
    d = DistancesAndCoordinatesLoader()
    d.import_distances(p.pdbs)
    d.import_coordinates(p.pdbs)
    """import info about redundant nucleotides"""
    r = RedundantNucleotidesLoader()
    r.import_redundant_nucleotides(p.pdbs)
    """import best chains and models"""
    b = BestChainsAndModelsLoader()
    b.import_best_chains_and_models(p.pdbs)

    exit()

    """cluster motifs"""
    c = ClusterMotifs()
    c.get_pdb_ids_for_clustering()
    c.get_loops_for_clustering(loop_type='IL')
    c.make_input_file_for_matlab()

#     c.cluster_loops(output_dir)

    # import motif atlas release into the database


    # import loop_searches, loop_positions and loop_searches_qa
    s = LoopSearchesLoader()
    s.load_loop_searches()
    s.load_loop_positions()
#     s.load_loop_search_qa_text_file('/Users/anton/FR3D/MM_extraNTs.txt')
#     s.load_loop_search_qa_text_file('/Users/anton/FR3D/MM_symmetrize.txt')


    # annotate all pdb files with these clusters

    # compute new non-redundant lists, import into the database



    logging.info('SUCCESSFUL UPDATE')

    # on failure: stop, email
    # log from matlab, log from python
    # log with dates, clear filenames


    m.send_report()


if __name__ == "__main__":
    main(sys.argv[1:])