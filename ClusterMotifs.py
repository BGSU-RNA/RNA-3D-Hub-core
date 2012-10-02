"""

Main entry point for motif clustering

"""

__author__ = 'Anton Petrov'

import os
import pdb
import sys
import getopt
import logging
import math
import time
from time import localtime, strftime
from subprocess import Popen, list2cmdline


from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import session, AllLoops, PdbBestChainsAndModels, NR_release, NR_pdb
from MotifLoader import MotifLoader


class ClusterMotifs(MotifAtlasBaseClass):
    """
    """
    def __init__(self, loop_type=None):
        MotifAtlasBaseClass.__init__(self)
        self.num_jobs = 4
        self.loop_type = loop_type
        self.pdb_ids = []
        self.loop_ids = []
        self.best_loops = [] # these loops will be clustered
        self.mlab_input_filename = os.path.join(os.getcwd(), 'loops.txt')
        self.__make_output_directory()

    def __make_output_directory(self):
        """make a directory for the release files"""
        self.output_dir = os.path.join( self.config['locations']['releases_dir'],
                                         self.loop_type + '_' + strftime("%Y%m%d_%H%M", localtime() ))
        if not os.path.exists(self.output_dir):
            os.mkdir( self.output_dir )
        logging.info('Files will be saved in %s' % self.output_dir)

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

    def get_loops_for_clustering(self):
        """get all loops from non-redundant pdbs"""
        loops = session.query(AllLoops). \
                        filter(AllLoops.pdb.in_(self.pdb_ids)). \
                        filter_by(type=self.loop_type). \
                        all()
        self.loop_ids = [loop.id for loop in loops]
        logging.info('Found %i loops' % len(loops) )

        """TODO: join with loop_qa table to filter out modified loops etc"""



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

    def parallel_exec_commands(self, cmds):
        """Exec commands in parallel in multiple process.
        Borrowed from:
        http://code.activestate.com/recipes/577376-simple-way-to-execute-multiple-process-in-parallel/
        """
        if not cmds:
            logging.critical('No commands to execute')
            return

        def done(p):
            return p.poll() is not None
        def success(p):
            return p.returncode == 0
        def fail():
            sys.exit(1)

        processes = []
        while True:
            while cmds and len(processes) < self.num_jobs:
                task = cmds.pop()
                list2cmdline(task)
                processes.append(Popen(task))

            for p in processes:
                if done(p):
                    if success(p):
                        processes.remove(p)
                    else:
                        fail()

            if not processes and not cmds:
                break
            else:
                time.sleep(0.05)

    def prepare_aAa_commands(self):
        """a list of matlab commands to run all-against-all searches in parallel
        """
        N = len(self.best_loops)
        interval = int(math.ceil( N / float(self.num_jobs)))
        logging.info('%i loops, will process in groups of %i' % (N, interval))
        commands = []
        mlab_app = '/Applications/MATLAB_R2007b/bin/matlab'
        mlab_params = '-nodisplay -nojvm -nodesktop -r '
        current_max = 0

        while current_max < N:
            aAa_command = "aAaSearches('%s', %i, %i)" % (self.mlab_input_filename,
                                                         current_max + 1,
                                                         current_max + interval)
            mlab_command = ';'.join(['"setup', aAa_command, 'quit"'])
            commands.append([mlab_app, mlab_params + mlab_command]) #"'setup; disp(rand(10)); quit();'"
            current_max += interval
        return commands

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

    M = ClusterMotifs(loop_type='IL')

    M.get_pdb_ids_for_clustering()
    M.get_loops_for_clustering()

#     M.manually_get_loops_for_clustering()

    M.make_input_file_for_matlab()

#     M.parallel_exec_commands( M.prepare_aAa_commands() )

    M.cluster_loops()
#
#     L = MotifLoader(motif_type='il')
#     L.import_data()


if __name__ == "__main__":
    main(sys.argv[1:])

