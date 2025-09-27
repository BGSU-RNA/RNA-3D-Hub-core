"""
The main entry point for motif clustering.
This used to run Matlab code, but now it simply runs Python code.

* Usage:

To run clustering on a single PDB file:
python ClusterMotifs.py 1FG0

To cluster all representative loops from the current non-redundant list:
python ClusterMotifs.py

"""

import os
# import math
# import time
import glob
from time import localtime, strftime
# from subprocess import Popen, list2cmdline

from pymotifs import core
from pymotifs import models as mod
from pymotifs.motif_atlas.compare_and_cluster import cluster_loops


class ClusterMotifs(core.Base):
    jobs = 4        # how many simultaneous Matlab jobs are going to be run
    script_prefix = 'aAa_script_'
    retries = 3

    def __init__(self, *args):
        super(ClusterMotifs, self).__init__(*args)
        self.fr3d_root = self.config['locations']['fr3d_root']
        self.mlab_input_filename = os.path.join(self.fr3d_root, 'loops.txt')

    def make_release_directory(self, loop_type, release_id):
        """
        Make a directory for the release files.
        The directory name will be based on the current time to avoid duplicatles.
        Directory name includes loop type and release number.

        Parameters
        ----------
        loop_type : str
            The loop type to create a directory for
        release_id : str
            The motif release we are creating

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

    # def make_input_file_for_matlab(self, loops, output_dir):
    #     # write for Matlab
    #     with open(self.mlab_input_filename, 'wb') as out:
    #         out.write(','.join(loops))

    #     self.logger.info('Saved loop_ids into %s' % self.mlab_input_filename)

    #     # store with the motif release data as well
    #     with open(os.path.join(output_dir, 'loops.txt'), 'wb') as out:
    #         out.write(','.join(loops))

    #     self.logger.info('Saved loop_ids into %s' % os.path.join(output_dir, 'loops.txt'))

    # def parallel_exec_commands(self, cmds):
    #     """Execute commands in parallel in multiple process.
    #     Adapted from:
    #     https://code.activestate.com/recipes/577376-simple-way-to-execute-multiple-process-in-parallel/

    #     If at least one of the parallel tasks fails and can't recover, the
    #     program will abort.
    #     """

    #     if not cmds:
    #         raise core.StageFailed("No commands to execute")

    #     processes = []
    #     tasks = {}
    #     retries_left = self.retries
    #     while True:
    #         while cmds and len(processes) < self.jobs:
    #             task = cmds.pop()
    #             list2cmdline(task)
    #             p = Popen(task)
    #             tasks[p.pid] = task  # associate task with a pid
    #             self.logger.info('Task %s has pid %i' % (task, p.pid))
    #             processes.append(p)

    #         for p in processes:
    #             if p.poll() is not None:
    #                 if p.returncode == 0:
    #                     self.logger.info('Parallel task succeeded')
    #                     processes.remove(p)
    #                 else:
    #                     retries_left -= 1
    #                     if not retries_left:
    #                         self.logger.critical('Parallel task failed')
    #                         self.logger.error(p.stdout)
    #                         self.logger.error(p.stderr)
    #                         raise matlab.MatlabFailed("Clustering failed")

    #                     self.logger.warning('Restarting Parallel task')
    #                     task = tasks[p.pid]  # retrieve the failed task
    #                     processes.remove(p)
    #                     p = Popen(task)  # new process with the same task
    #                     tasks[p.pid] = task  # save the new task's pid
    #                     processes.append(p)

    #         if not processes and not cmds:
    #             break
    #         else:
    #             time.sleep(0.05)

    # def prepare_aAa_commands(self, loops, enforceSize=True):
    #     """Creates a list of matlab commands to run all-against-all searches
    #     in parallel. To avoid matlab hanging at the command prompt in case
    #     of errors in the matlab code, the script must be written out to a
    #     file, and then this file is launched wrapped up in a try/catch
    #     statement.
    #     """

    #     N = len(loops)
    #     interval = int(math.ceil(N / float(self.jobs)))
    #     self.logger.info('%i loops, will process in groups of %i' %
    #                      (N, interval))
    #     commands = []
    #     mlab_params = ' -nodisplay -nojvm -r '
    #     current_max = 0

    #     base = self.config['locations']['base']
    #     i = 1
    #     while current_max < N:
    #         # prepare matlab code
    #         mlab_command = SCRIPT.format(base=base,
    #                                      fr3d=self.fr3d_root,
    #                                      input_file=self.mlab_input_filename,
    #                                      start=current_max + 1,
    #                                      stop=current_max + interval,
    #                                      enforceSize=int(enforceSize))

    #         script_name = '%s%i.m' % (self.script_prefix, i)
    #         i += 1

    #         # save matlab code to a temporary matlab script
    #         script_path = os.path.join(self.fr3d_root, script_name)
    #         with open(script_path, 'wb') as out:
    #             out.write(mlab_command)

    #         # prepare bash command to run the matlab script
    #         ext = script_name[:-2]
    #         try_catch = 'try %s(), catch, exit(1), end, exit(0)' % ext
    #         cd_command = "cd {fr3d}".format(fr3d=self.fr3d_root)
    #         bash_command = '"' + ';'.join([cd_command, try_catch]) + '"'
    #         commands.append([self.config['locations']['mlab_app'],
    #                         mlab_params + bash_command])
    #         current_max += interval

    #     for cmd in commands:
    #         self.logger.info(cmd[1])

    #     return commands

    def _clean_up(self):
        """Remove temporary file with loop ids and temporary .m scripts
        """

        if os.path.exists(self.mlab_input_filename):
            os.remove(self.mlab_input_filename)

        scripts = os.path.join(self.fr3d_root, self.script_prefix + '*.m')

        for filename in glob.glob(scripts):
            os.remove(filename)

    def __call__(self, loop_type, loop_position_to_border_unit_id, release_id, molecule_type):
        """
        Cluster motifs of the given type for the given pdb files.
        This will get all
        valid loops from the best chains and models and then run a series of
        parallel matlab jobs to cluster them.

        loops is really loop_position_border_units, a complicated dictionary

        :loop_type: The type of loops to cluster.
        :pdbs: The pdb files to get loops from.
        :returns: Nothing.
        """

        # not sure why this is here
        import shutil
        def copytree(src, dst):
            # List contents of source directory
            for item in os.listdir(src):
                # Construct full path to the item
                source_item_path = os.path.join(src, item)
                destination_item_path = os.path.join(dst, item)

                # If item is a directory, copy its contents (recursive call)
                if os.path.isdir(source_item_path):
                    # If a directory with the same name doesn't exist in the destination, create a new one
                    if not os.path.exists(destination_item_path):
                        shutil.copytree(source_item_path, destination_item_path)
                    else:
                        # If the directory already exists, recursively copy into it
                        copytree(source_item_path, destination_item_path)
                else:
                    # It's a file, copy it directly
                    shutil.copy2(source_item_path, destination_item_path)  # copy2 preserves metadata
                print("Copied %s to %s" % (source_item_path, destination_item_path))


        if not loop_position_to_border_unit_id:
            raise ValueError("Must give loops to cluster")

        self._clean_up()

        output_dir = self.make_release_directory(loop_type, release_id)

        # retrieve the table of loop annotations and write those to the release directory
        loop_annotations = []
        with self.session() as session:
            query = session.query(mod.LoopAnnotations)

            for row in query:
                text = "%s\t%s\t%s\t%s" % (row.loop_id, row.annotation_1, row.annotation_2, row.author)
                loop_annotations.append(text)
        with open(os.path.join(output_dir, 'loop_annotations.txt'), 'w') as out:
            out.write('\n'.join(loop_annotations))

        # cluster the current loop_type
        print('Running cluster_loops method of compare_and_cluster.py')
        self.logger.info('Running cluster_loops method of compare_and_cluster.py')

        # must each motif group have at least one loop solved by x-ray?
        # get the list of PDB ids that are solved by x-ray
        x_ray_pdb_list = []
        with self.session() as session:
            # get PDB ids with experimental technique including 'X-RAY DIFFRACTION'
            query = session.query(mod.PdbInfo.pdb_id).\
                filter(mod.PdbInfo.experimental_technique.like('%X-RAY DIFFRACTION%'))

            x_ray_pdb_list = [row.pdb_id for row in query]

        # self.logger.info('X-ray PDB list: %s' % x_ray_pdb_list)

        cluster_loops(loop_position_to_border_unit_id, output_dir, molecule_type, x_ray_pdb_list)

        self.logger.info('Successful clustering of %s' % loop_type)

        self._clean_up()

        self.logger.info('Cleaned up')

        return output_dir
