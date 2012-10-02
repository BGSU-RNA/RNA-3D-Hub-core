"""

Program for running all-against-all searches in parallel.

"""

import sys
import os
import time
import logging
from subprocess import Popen

from MotifAtlasBaseClass import MotifAtlasBaseClass


class ParallelAllAgainstAll(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        """
        """
        self.num_jobs = 4

    def exec_commands(self, cmds):
        """Exec commands in parallel in multiple process"""
        if not cmds: return # empty list

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


def main(argv):
    """
    """

    logging.basicConfig(level=logging.DEBUG)

    loop_ids = []

    P = ParallelAllAgainstAll(loop_ids)

    commands = [
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'setup;aAaSearches('filename', 1, 100);'"],
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'disp(rand(10));quit();'"],
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'disp(rand(10));quit();'"],
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'disp(rand(10));quit();'"],
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'disp(rand(10));quit();'"],
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'disp(rand(10));quit();'"],
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'disp(rand(10));quit();'"],
        ["/Applications/MATLAB_R2007b/bin/matlab", "-nodisplay -nojvm -nodesktop -r 'disp(rand(10));quit();'"]
    ]
    P.exec_commands(commands)



if __name__ == "__main__":
    main(sys.argv[1:])

