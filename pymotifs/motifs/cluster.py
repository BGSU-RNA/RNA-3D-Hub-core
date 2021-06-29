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

If matlab fails, python tries to restart the same process several times. This
is done based on the observation that sometimes a FR3D search can crash, but
will run smoothly the next time round without any intervention.
"""

import os
import math
import time
import glob
from time import localtime, strftime
from subprocess import Popen, list2cmdline

from pymotifs import core
from pymotifs.utils import matlab

SCRIPT = """
cd '{base}'
setup()
cd '{fr3d}'
aAaSearches('{input_file}', {start}, {stop}, {enforceSize})
"""


# It may be useful to log to the log file using:
# https://codereview.stackexchange.com/questions/6567/redirecting-subprocesses-output-stdout-and-stderr-to-the-logging-module

class ClusterMotifs(core.Base):
    jobs = 4
    script_prefix = 'aAa_script_'
    retries = 3

    def __init__(self, *args):
        super(ClusterMotifs, self).__init__(*args)
        self.fr3d_root = self.config['locations']['fr3d_root']
        self.mlab_input_filename = os.path.join(self.fr3d_root, 'loops.txt')

    def make_release_directory(self, loop_type, release_id):
        """Make a directory for the release files.
        The directory name will be based on the current time to avoid duplicatles.
        Directory name includes loop type and release number.

        Parameters
        ----------
        loop_type : str
            The loop type to create a directory for.
        release_id : str
            The representative release we are working with.

        Returns
        -------
        directory : str
            Full path to the created directory.
        """

        release_dir = self.config['locations']['releases_dir']
        time_stamp = strftime("%Y-%m-%d_%H:%M", localtime())
        output_dir = os.path.join(release_dir, loop_type + '_' + release_id + '_' + time_stamp)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        else:
            self.logger.warning("Motif directory already exists")
        self.logger.info('Files will be saved in %s' % output_dir)
        return output_dir

    def make_input_file_for_matlab(self, loops):
        with open(self.mlab_input_filename, 'wb') as out:
            out.write(','.join(loops))

        self.logger.info('Saved loop_ids into %s' % self.mlab_input_filename)

    def parallel_exec_commands(self, cmds):
        """Execute commands in parallel in multiple process.
        Adapted from:
        http://code.activestate.com/recipes/577376-simple-way-to-execute-multiple-process-in-parallel/

        If at least one of the parallel tasks fails and can't recover, the
        program will abort.
        """

        if not cmds:
            raise core.StageFailed("No commands to execute")

        processes = []
        tasks = {}
        retries_left = self.retries
        while True:
            while cmds and len(processes) < self.jobs:
                task = cmds.pop()
                list2cmdline(task)
                p = Popen(task)
                tasks[p.pid] = task  # associate task with a pid
                self.logger.info('Task %s has pid %i' % (task, p.pid))
                processes.append(p)

            for p in processes:
                if p.poll() is not None:
                    if p.returncode == 0:
                        self.logger.info('Parallel task succeeded')
                        processes.remove(p)
                    else:
                        retries_left -= 1
                        if not retries_left:
                            self.logger.critical('Parallel task failed')
                            self.logger.error(p.stdout)
                            self.logger.error(p.stderr)
                            raise matlab.MatlabFailed("Clustering failed")

                        self.logger.warning('Restarting Parallel task')
                        task = tasks[p.pid]  # retrieve the failed task
                        processes.remove(p)
                        p = Popen(task)  # new process with the same task
                        tasks[p.pid] = task  # save the new task's pid
                        processes.append(p)

            if not processes and not cmds:
                break
            else:
                time.sleep(0.05)

    def prepare_aAa_commands(self, loops, enforceSize=True):
        """Creates a list of matlab commands to run all-against-all searches
        in parallel. To avoid matlab hanging at the command prompt in case
        of errors in the matlab code, the script must be written out to a
        file, and then this file is launched wrapped up in a try/catch
        statement.
        """

        N = len(loops)
        interval = int(math.ceil(N / float(self.jobs)))
        self.logger.info('%i loops, will process in groups of %i' %
                         (N, interval))
        commands = []
        mlab_params = ' -nodisplay -nojvm -r '
        current_max = 0

        base = self.config['locations']['base']
        i = 1
        while current_max < N:
            # prepare matlab code
            mlab_command = SCRIPT.format(base=base,
                                         fr3d=self.fr3d_root,
                                         input_file=self.mlab_input_filename,
                                         start=current_max + 1,
                                         stop=current_max + interval,
                                         enforceSize=int(enforceSize))

            script_name = '%s%i.m' % (self.script_prefix, i)
            i += 1

            # save matlab code to a temporary matlab script
            script_path = os.path.join(self.fr3d_root, script_name)
            with open(script_path, 'wb') as out:
                out.write(mlab_command)

            # prepare bash command to run the matlab script
            ext = script_name[:-2]
            try_catch = 'try %s(), catch, exit(1), end, exit(0)' % ext
            cd_command = "cd {fr3d}".format(fr3d=self.fr3d_root)
            bash_command = '"' + ';'.join([cd_command, try_catch]) + '"'
            commands.append([self.config['locations']['mlab_app'],
                            mlab_params + bash_command])
            current_max += interval

        for cmd in commands:
            self.logger.info(cmd[1])

        return commands

    def _clean_up(self):
        """Remove temporary file with loop ids and temporary .m scripts
        """

        if os.path.exists(self.mlab_input_filename):
            os.remove(self.mlab_input_filename)

        scripts = os.path.join(self.fr3d_root, self.script_prefix + '*.m')

        for filename in glob.glob(scripts):
            os.remove(filename)

    def __call__(self, loop_type, loops, release_id):
        """Launch the main matlab motif clustering pipeline. This will cluster
        motifs of the given type for the given pdb files. This will get all
        valid loops from the best chains and models and then run a series of
        parallel matlab jobs to cluster them.

        :loop_type: The type of loops to cluster.
        :pdbs: The pdb files to get loops from.
        :returns: Nothing.
        """

        if not loops:
            raise ValueError("Must give loops to cluster")

        self._clean_up()

        output_dir = self.make_release_directory(loop_type, release_id)
        self.make_input_file_for_matlab(loops)
        self.parallel_exec_commands(self.prepare_aAa_commands(loops))

        mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
        [status, err_msg] = \
            mlab.MotifAtlasPipeline(self.mlab_input_filename, output_dir, nout=2)

        if err_msg:
            raise matlab.MatlabFailed(err_msg)

        self._clean_up()
        self.logger.info('Successful clustering of %s' % loop_type)

        return output_dir
