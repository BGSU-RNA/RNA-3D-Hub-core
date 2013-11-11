"""

The main entry point for motif clustering.


* Usage:

To run clustering on a single PDB file:
python ClusterMotifs.py 1FG0

To cluster all representative loops from the current non-redundant list:
python ClusterMotifs.py


* Notes on parallelization:

The program will split all-against-all searches into an approximately equal
batch jobs and submit them to several Matlab processes running in parallel.

In most other programs supporting RNA 3D Hub, Matlab is called via mlabwrap. In
this program, however, matlab is called directly using subprocesses. This is
done in order to parallelize and speed up all-against-all searches, but it is
inherently less reliable.

Although other parts of the RNA 3D Hub pipeline also yield themselves to
parallelization, it might be more appropriate to use it in this context because
motif clustering will be run less often than the rest of the pipeline, and the
potential errors will have smaller impact.

If matlab fails, python tries to restart the same process several times. This is
done based on the observation that sometimes a FR3D search can crash, but will
run smoothly the next time round without any intervention.

"""

__author__ = 'Anton Petrov'

import os
import pdb
import sys
import getopt
import logging
import shutil
import math
import time
import glob
import pdb
import datetime
from time import localtime, strftime
from subprocess import Popen, list2cmdline


from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import session, AllLoops, PdbBestChainsAndModels, NR_release, \
                   NR_pdb, LoopRelease, LoopQA, Release
from MotifLoader import MotifLoader


class ClusterMotifs(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.success    = False
        self.num_jobs   = 4
        self.pdb_ids    = []
        self.loop_ids   = []
        self.best_loops = [] # loops to be clustered
        self.fr3d_root  = self.config['locations']['fr3d_root']
        self.retries_left  = 3
        self.script_prefix = 'aAa_script_'
        self.mlab_input_filename = os.path.join(self.fr3d_root, 'loops.txt')

    def set_loop_type(self, loop_type):
        """
        """
        self.loop_type = loop_type

    def is_four_weeks_since_last_update(self):
        """
            Get the latest release of the correct loop type and see if the
            specified amount of time has elapsed since it was produced.
            Return true or false.
        """
        previous = session.query(Release).filter(Release.type==self.loop_type).\
                                          order_by(Release.date.desc()).\
                                          first()
        logging.info('Last %s release occurred on %s'
                          % (self.loop_type,
                             datetime.datetime.strftime(previous.date, "%Y-%m-%d %H:%M")))
        # 27 days is after 3 weekly update cycles
        if (previous.date + datetime.timedelta(days=27)) > datetime.datetime.now():
            logging.info('Next %s release will occur on %s, now skipping.' \
                          % (self.loop_type,
                             datetime.datetime.strftime(previous.date + datetime.timedelta(weeks=4), "%Y-%m-%d %H:%M")))
            return False
        else:
            logging.info('Time for new %s clustering' % self.loop_type)
            return True

    def make_release_directory(self):
        """make a directory for the release files"""
        self.output_dir = os.path.join( self.config['locations']['releases_dir'],
                                        self.loop_type + '_' + strftime("%Y%m%d_%H%M", localtime() ))
        if not os.path.exists(self.output_dir):
            os.makedirs( self.output_dir )
        logging.info('Files will be saved in %s' % self.output_dir)

    def _remove_release_directory(self):
        """Especially useful for tests"""
        shutil.rmtree(self.output_dir)

    def get_pdb_ids_for_clustering(self, nr_release_id=None):
        """
            Use the specified NR release or the latest one by default.
        """
        if nr_release_id:
            latest_nr_release = session.query(NR_release). \
                                        filter(NR_release.id==nr_release_id). \
                                        first()
        else:
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
        """ideally this should be achieved with joins, but it requires adjusting
        the models"""
        """get latest loop release"""
        latest_loop_release = session.query(LoopRelease).\
                                      order_by(LoopRelease.date.desc()).\
                                      first()
        """get all valid loops"""
        valid_ids = [loop.id for loop in session.query(LoopQA).\
                              filter(LoopQA.status==1).\
                              filter(LoopQA.release_id==latest_loop_release.id).\
                              filter(LoopQA.id.like(self.loop_type + '%')).\
                              all()]
        """get all loops from non-redundant pdbs"""
        loops = session.query(AllLoops). \
                        filter(AllLoops.pdb.in_(self.pdb_ids)). \
                        filter_by(type=self.loop_type). \
                        all()
        logging.info('Found %i loops' % len(loops))
        """filter out invalid loops"""
        loops = [loop for loop in loops if loop.id in valid_ids]
        self.loop_ids = [loop.id for loop in loops if loop.id in valid_ids]
        logging.info('Kept %i valid loops' % len(loops))
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
        """remove manually blacklisted loops. In future use database"""
        blacklist = ['HL_3ICQ_004', 'HL_3V2F_005', 'HL_1Y0Q_002', 'HL_2IL9_002', 'HL_2IL9_005', 'HL_3V2F_065']
        for bad_loop in blacklist:
            if bad_loop in self.best_loops:
                self.best_loops.remove(bad_loop)
                logging.info('Removed blacklisted loop %s' % bad_loop)
        logging.info('Selected %i loops', len(self.best_loops))

    def make_input_file_for_matlab(self):
        """
        """
        f = open(self.mlab_input_filename, 'w')
        f.write(','.join(self.best_loops))
        f.close()
        logging.info('Saved loop_ids in file %s' % self.mlab_input_filename)

    def parallel_exec_commands(self, cmds):
        """
            Exec commands in parallel in multiple process.
            Adapted from:
            http://code.activestate.com/recipes/577376-simple-way-to-execute-multiple-process-in-parallel/

            It at least one of the parallel tasks fails and can't recover, the
            program will abort.
        """
        if not cmds:
            logging.critical('No commands to execute')
            self._crash()

        def done(p):
            return p.poll() is not None
        def success(p):
            return p.returncode == 0

        processes = []
        tasks = {}
        while True:
            while cmds and len(processes) < self.num_jobs:
                task = cmds.pop()
                list2cmdline(task)
                p = Popen(task)
                tasks[p.pid] = task # associate task with a pid
                logging.info('%i %s' % (p.pid, task))
                processes.append(p)

            for p in processes:
                if done(p):
                    if success(p):
                        logging.info('Parallel task completed successfully')
                        processes.remove(p)
                    else:
                        self.retries_left -= 1
                        if self.retries_left == 0:
                            logging.critical('Parallel task failed')
                            self._crash()
                        else:
                            logging.warning('Parallel task will be restarted')
                            task = tasks[p.pid] # retrieve the failed task
                            processes.remove(p)
                            p = Popen(task) # new process with the same task
                            tasks[p.pid] = task # save the new task's pid
                            processes.append(p)

            if not processes and not cmds:
                break
            else:
                time.sleep(0.05)

    def prepare_aAa_commands(self):
        """
            Creates a list of matlab commands to run all-against-all searches
            in parallel. To avoid matlab hanging at the command prompt in case
            of errors in the matlab code, the script must be written out to a
            file, and then this file is launched wrapped up in a try/catch
            statement.
        """
        N = len(self.best_loops)
        interval = int(math.ceil( N / float(self.num_jobs)))
        logging.info('%i loops, will process in groups of %i' % (N, interval))
        commands = []
        mlab_params = ' -nodisplay -nojvm -r '
        current_max = 0

        os.chdir(self.fr3d_root)

        i = 1
        while current_max < N:
            # prepare matlab code
            cd_command = "cd '%s'" % self.fr3d_root
            aAa_command = "aAaSearches('%s', %i, %i)" % (self.mlab_input_filename,
                                                         current_max + 1,
                                                         current_max + interval)
            mlab_command = ';'.join([cd_command, 'setup', aAa_command])
            script_name = '%s%i.m' % (self.script_prefix, i)
            i += 1
            # save matlab code to a temporary m script
            script_path = os.path.join(self.fr3d_root, script_name)
            f = open(script_path, 'w')
            f.write(mlab_command)
            f.close()
            # prepare bash command to run the matlab script
            try_catch  = 'try %s(), catch, exit(1), end, exit(0)' % script_name[:-2]
            bash_command = '"' + ';'.join([cd_command, try_catch]) + '"'
            commands.append([self.config['locations']['mlab_app'],
                            mlab_params + bash_command])
            current_max += interval
        [logging.info(x[1]) for x in commands]
        return commands

    def _clean_up(self):
        """
            remove temporary file with loop ids and temporary .m scripts
        """
        os.remove(self.mlab_input_filename)
        temp_scripts = os.path.join(self.fr3d_root,
                                    self.script_prefix + '*.m')
        [os.remove(x) for x in glob.glob( temp_scripts )]

    def cluster_loops(self):
        """
            Launch the main matlab motif clustering pipeline.
        """
        try:
            self._setup_matlab()
            [status, err_msg] = self.mlab.MotifAtlasPipeline(self.mlab_input_filename,
                                                             self.output_dir,
                                                             nout=2)
            self._clean_up()
            if err_msg == '':
                self.success = True
                logging.info('Successful clustering')
            else:
                logging.critical(err_msg)
        except:
            e = sys.exc_info()[1]
            self._crash(e)

    def _manually_get_loops_for_clustering(self, pdb_id):
        """
            For testing purposes.
        """
        self.best_loops = [loop.id for loop in session.query(AllLoops).
                                                       filter(AllLoops.pdb==pdb_id).
                                                       filter(AllLoops.type=='il').
                                                       all()]
        logging.info('Selected %s loops from %s' % (len(self.best_loops), pdb_id) )


def main(argv):
    """
    """

    M = ClusterMotifs()
    M.start_logging()
    M.set_loop_type('HL')
    M.make_release_directory()

    if len(argv) > 0:
        M._manually_get_loops_for_clustering(argv[0])
    else:
        M.get_pdb_ids_for_clustering()
        M.get_loops_for_clustering()

    M.make_input_file_for_matlab()

    M.parallel_exec_commands( M.prepare_aAa_commands() )

    M.cluster_loops()

    M.set_email_subject('Successful clustering')
    M.send_report()


if __name__ == "__main__":
    main(sys.argv[1:])
