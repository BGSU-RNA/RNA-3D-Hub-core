"""

Main entry point for launching a Motif Atlas update.

Does not yet include NR list update and motif clustering.

Examples:

* to launch the pipeline:
    python MotifAtlas.py

* to run just on one test pdb file:
    python MotifAtlas.py test

"""

__author__ = 'Anton Petrov'

import sys
import logging
import pdb

from DistancesAndCoordinatesLoader import DistancesAndCoordinatesLoader
from PdbInfoLoader import PdbInfoLoader
from LoopExtractor import LoopExtractor
from LoopQualityChecker import LoopQualityChecker
from PairwiseInteractionsLoader import PairwiseInteractionsLoader
from RedundantNucleotidesLoader import RedundantNucleotidesLoader
from BestChainsAndModelsLoader import BestChainsAndModelsLoader
from MotifAtlasBaseClass import MotifAtlasBaseClass
from ClusterMotifs import ClusterMotifs
from LoopSearchesLoader import LoopSearchesLoader
from CacheManager import CacheManager
from PdbFileExporter import PdbFileExporter


def main(argv):

    try:
        """set up logging"""
        m = MotifAtlasBaseClass()
        m.start_logging()
        logging.info('Initializing update')

        """get new pdb files, import descriptions into the database"""
        p = PdbInfoLoader()

        if argv and argv[0] == 'test':
            p.pdbs = ['1FG0']
        else:
            p.get_all_rna_pdbs()
            p.update_rna_containing_pdbs()
            p.check_obsolete_structures()

        """extract all loops and import into the database"""
        e = LoopExtractor()
        e.extract_and_import_loops(p.pdbs)

        """do loop QA, import into the database.
        This will always create a new loop release."""
        q = LoopQualityChecker()
        q.check_loop_quality(p.pdbs)

         """import pairwise interactions annotated by FR3D"""
        i = PairwiseInteractionsLoader()
        i.import_interactions(p.pdbs)

        """export pairwise interactions to a compressed file for NDB"""
        f = PdbFileExporter()
        f.export_interactions(m.config['locations']['interactions_gz'])

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

        """update cache"""
        c = CacheManager('rna3dhub')
        c.update_pdb_cache()
        c.update_nrlist_cache()
        c.update_motif_cache()
        c.update_loop_cache()

        logging.info('SUCCESSFUL UPDATE')

        m.send_report()

        """cluster motifs"""
        #     c = ClusterMotifs(loop_type='IL')
        #     c.get_pdb_ids_for_clustering()
        #     print c.pdb_ids
        #     c.pdb_ids = ['1S72']
        #     c.get_loops_for_clustering()
        #     c.make_input_file_for_matlab()
        #     c.cluster_loops()
        # import motif atlas release into the database
        # import loop_searches, loop_positions and loop_searches_qa
        #     s = LoopSearchesLoader()
        #     s.load_loop_searches()
        #     s.load_loop_positions()
        #     s.load_loop_search_qa_text_file('FR3D/MM_extraNTs.txt')
        #     s.load_loop_search_qa_text_file('FR3D/MM_symmetrize.txt')
        # annotate all pdb files with these clusters
        # compute new non-redundant lists, import into the database

    except:
        try:
            logging.critical('Update FAILED')
            m.send_report()
        except:
            pass
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])