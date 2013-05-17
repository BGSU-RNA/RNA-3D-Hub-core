"""

Main entry point for RNA 3D Hub updates.

Examples:

* to launch the pipeline:
    python MotifAtlas.py

"""

__author__ = 'Anton Petrov'

import sys
import os
import logging
import traceback
import pdb


from BestChainsAndModelsLoader import BestChainsAndModelsLoader
from CacheManager import CacheManager
from ClusterMotifs import ClusterMotifs
from DistancesAndCoordinatesLoader import DistancesAndCoordinatesLoader
from LoopExtractor import LoopExtractor
from LoopQualityChecker import LoopQualityChecker
from LoopSearchesLoader import LoopSearchesLoader
from MotifAtlasBaseClass import MotifAtlasBaseClass
from MotifLoader import MotifLoader
from NrListLoader import NrListLoader
from PairwiseInteractionsLoader import PairwiseInteractionsLoader
from PdbDownloader import PdbDownloader
from PdbFileExporter import PdbFileExporter
from PdbInfoLoader import PdbInfoLoader
from RedundantNucleotidesLoader import RedundantNucleotidesLoader
from unit_ordering_loader import UnitOrderingLoader
from UnitIdLoader import UnitIdLoader


def update_unit_ordering(pdb_ids):
    loader = UnitOrderingLoader()
    try:
        loader.import_ordering(pdb_ids)
    except:
        logging.error(traceback.format_exc(sys.exc_info()))
        logging.error('Could not compute ordering')


def cluster_motifs(motif_type, nr_release_id=None):
    """
        cluster motifs
    """
    try:
        c = ClusterMotifs()
        c.set_loop_type(motif_type)
        if not c.is_four_weeks_since_last_update():
            return
        c.make_release_directory()
        c.get_pdb_ids_for_clustering(nr_release_id)
        c.get_loops_for_clustering()
        c.make_input_file_for_matlab()
        c.parallel_exec_commands( c.prepare_aAa_commands() )
        c.cluster_loops()

        l = LoopSearchesLoader()
        l.load_loop_search_qa_text_file(os.path.join(c.output_dir, 'MM_extraNTs.txt'))
        l.load_loop_search_qa_text_file(os.path.join(c.output_dir, 'MM_symmetrize.txt'))

        logging.info('%s loops successfully clustered' % motif_type)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('%s clustering failed' % motif_type)


def import_motifs(motif_type):
    """
        import loop clusters into the database
    """
    try:
        m = MotifLoader(motif_type=motif_type)
        m.import_data()
        logging.info('%s loops imported' % motif_type)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('%s loop import failed' % motif_type)

def update_pdb_info():
    """
        get new pdb files, import descriptions into the database.
        Return an empty list if something goes wrong, then all following
        steps that rely on pdb_ids will be skipped.
    """
    try:
        p = PdbInfoLoader()
        p.get_all_rna_pdbs()
        pdb_ids = p.pdbs
        p.update_pdb_info()
        p.check_obsolete_structures()
#         p.update_organism_names()
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        pdb_ids = []

    return pdb_ids

def update_loops(pdb_ids):
    """
    """
    try:
        """extract all loops and import into the database"""
        e = LoopExtractor()
        e.extract_and_import_loops(pdb_ids)
        """do loop qa, create a new loop release. It only makes sense
        to do this if loop extraction worked"""
        q = LoopQualityChecker()
        q.check_loop_quality(pdb_ids)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Loop extraction or loop QA failed')

def update_pairwise_annotations(pdb_ids):
    """
    """
    try:
        """import pairwise interactions annotated by FR3D"""
        i = PairwiseInteractionsLoader()
        i.import_interactions(pdb_ids)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Pairwise interactions import failed')

def update_unit_ids(pdb_ids):
    """
        create new-style ids, must be done before exporting any data to files.
    """
    try:
        u = UnitIdLoader()
        u.import_unit_ids(pdb_ids)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Unit id import failed')

def update_cache():
    """
        update CodeIgniter cache
    """
    try:
        c = CacheManager()
    except:
        logging.warning("Skipping CodeIgniter cache update")
        return

    try:
        c.update_pdb_cache()
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Pdb cache update failed')

    try:
        c.update_nrlist_cache()
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Nrlist cache update failed')

    try:
        c.update_motif_cache()
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Motif cache update failed')

    try:
        c.update_loop_cache()
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Loop cache update failed')

def export_data(pdb_ids):
    """
        export data to static compressed files.
    """
    try:
        """export pairwise interactions to a compressed file for NDB"""
        f = PdbFileExporter()
        f.export_interactions(f.config['locations']['interactions_gz'], pdb_ids)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Pairwise interactions export failed')

def update_coordinates(pdb_ids):
    """
        import coordinates and distances into the database
    """
    try:
        d = DistancesAndCoordinatesLoader()
        d.import_distances(pdb_ids)
        d.import_coordinates(pdb_ids)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Distances or coordinates import failed')

def update_redundant_nucleotides(pdb_ids):
    """
        import info about redundant nucleotides
    """
    try:
        r = RedundantNucleotidesLoader()
        r.import_redundant_nucleotides(pdb_ids)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Redundant nt import failed')

def update_best_chains_and_models(pdb_ids):
    """
        import best chains and models
    """
    try:
        b = BestChainsAndModelsLoader()
        b.import_best_chains_and_models(pdb_ids)
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Redundant nt import failed')

def update_loop_annotations():
    """
        Must follow loop clustering because it imports all-against-all search
        data.
    """
    try:
        """import loop-level annotations"""
        LoopLoader = LoopSearchesLoader()
        LoopLoader.load_loop_positions()
        """import all-against-all search results from IL and HL clustering"""
        LoopLoader.load_loop_searches()
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Hairpin loop import failed')

def update_nrlists():
    """
        import non-redundant lists of RNA 3D structures.
    """
    try:
        """update Matlab files"""
        n = NrListLoader()
        """create file with resolution and source organism info"""
#         n.make_report_file()
#         n.update_nrlists()
        """create text files with equivalence classes"""
#         n.generate_output_files()
        """clean up"""
#         n.archive_files()
        """import into the database"""
        n.import_data()
        """create old style output using Matlab"""
#         n.make_old_style_html_tables()
#         n.update_old_website()
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('NR list update failed')

def download_files(pdb_ids):
    """
        download pdb and cif files from PDB. This has to be done independently
        from Matlab because there were connectivity problems between PDB and
        Matlab R2007b in the past.
    """
    try:
        d = PdbDownloader()
        d.pdbs = pdb_ids
        d.set_locations([d.config['locations']['cif_dir']])
        d.download_files()
        logging.info('Pdb and cif files successfully downloaded')
    except:
        logging.warning(traceback.format_exc(sys.exc_info()))
        logging.warning('Downloading pdb and cif files failed')


def main(argv):
    """
        RNA 3D Hub update entry point.
    """

    try:
        m = MotifAtlasBaseClass()
        m.start_logging()

#         pdb_ids = update_pdb_info()

#         download_files(pdb_ids)

        # must follow update_pdb_info
        update_nrlists()

        pdb.set_trace()

        update_loops(pdb_ids)

        update_pairwise_annotations(pdb_ids)

        update_unit_ids(pdb_ids)

        update_unit_ordering(pdb_ids)

        # must follow unit id updates
        export_data(pdb_ids)

        update_coordinates(pdb_ids)

        update_redundant_nucleotides(pdb_ids)

        update_best_chains_and_models(pdb_ids)

        # must follow best chain and model update
        cluster_motifs('IL')
        import_motifs('IL')

        cluster_motifs('HL')
        import_motifs('HL')

        # must follow motif clustering
        update_loop_annotations()

        # TODO annotate all pdb files with motifs

        update_cache()

        logging.info('Update completed')
        m.send_report()

    except:
        try:
            logging.critical('Update failed')
            logging.critical(traceback.format_exc(sys.exc_info()))
            m.set_email_subject('RNA 3D Hub update failed')
            m.send_report()
        except:
            pass


if __name__ == "__main__":
    main(sys.argv[1:])
