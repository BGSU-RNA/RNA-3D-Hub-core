"""

Main entry point for clustering motifs.

"""

__author__ = 'Anton Petrov'

import os
import pdb
import sys
import getopt
import logging
from time import localtime, strftime

from MotifAtlasBaseClass import MotifAtlasBaseClass
from MLSqlAlchemyClasses import session, AllLoops, PdbBestChainsAndModels
from NRSqlAlchemyClasses import NR_release, NR_pdb
from MLMotifLoader import Loader


class ClusterMotifs(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.pdb_ids = []
        self.loop_ids = []
        self.best_loops = [] # these loops will be clustered
        self.mlab_input_filename = os.path.join(os.getcwd(), 'loops.txt')
        """make a directory for the release files"""
        self.output_dir = os.path.join( self.config['locations']['releases_dir'],
                                         strftime("%Y%m%d_%H%M", localtime() ))
        print self.output_dir
        os.mkdir( self.output_dir )


    def get_pdb_ids_for_clustering(self):
        """get latest nr release"""
        latest_nr_release = session.query(NR_release). \
                                    order_by(NR_release.date.desc()). \
                                    first()
        logging.info('Will use NR release %s' % latest_nr_release.id)
        """get all pdbs from that nr release"""
        pdbs = session.query(NR_pdb). \
                       filter_by(release_id=latest_nr_release.id). \
                       filter_by(rep=1). \
                       filter(NR_pdb.class_id.like('NR_4%')). \
                       all()
        logging.info('Found %i NR pdbs' % len(pdbs))
        self.pdb_ids = [x.id for x in pdbs]

    def get_loops_for_clustering(self, loop_type):
        """get all loops from non-redundant pdbs"""
        loops = session.query(AllLoops). \
                        filter(AllLoops.pdb.in_(self.pdb_ids)). \
                        filter_by(type=loop_type). \
                        all()
        self.loop_ids = [loop.id for loop in loops]
        logging.info('Found %i loops' % len(loops) )
        """get info about best chains"""
        best_chains = dict()
        for x in session.query(PdbBestChainsAndModels).all():
            best_chains[x.pdb_id] = x.best_chains
        """keep only loops from best chains based on their nt_ids"""
        for loop in loops:
            chains = ''
            for nt_id in loop.nt_ids.split(','):
                chains += nt_id.split('_')[3]
            if loop.pdb in best_chains and list(set(chains) & set(best_chains[loop.pdb])):
                self.best_loops.append(loop.id)
                logging.info('Loop %s from chains %s belongs to best chains %s' \
                % (loop.id, chains, best_chains[loop.pdb] ))
        logging.info('Selected %i loops', len(self.best_loops))

    def make_input_file_for_matlab(self):
        """
        """
        f = open(self.mlab_input_filename, 'w')
        f.write(','.join(self.best_loops))
        f.close()
        logging.info('Saved loop_ids in file %s' % self.mlab_input_filename)

    def cluster_loops(self):
        """
        """
        try:
            MotifAtlasBaseClass._setup_matlab(self)
            [status, err_msg] = self.mlab.MotifAtlasPipeline(self.mlab_input_filename,
                                                             self.output_dir,
                                                             nout=2)
            os.remove(self.mlab_input_filename)
        except:
            e = sys.exc_info()[1]
            MotifAtlasBaseClass._crash(self,e)

    def manually_get_loops_for_clustering(self):
        """
        """
        self.best_loops = [loop.id for loop in session.query(AllLoops).
                                                       filter(AllLoops.pdb=='1S72').
                                                       filter(AllLoops.type=='il').
                                                       all()]


def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)

    M = ClusterMotifs()

#     M.get_pdb_ids_for_clustering()
#     M.get_loops_for_clustering(loop_type='IL')

    M.manually_get_loops_for_clustering()

    M.make_input_file_for_matlab()

    M.cluster_loops()
#
#     L = Loader(motif_type='il')
#     L.import_data()


if __name__ == "__main__":
    main(sys.argv[1:])

