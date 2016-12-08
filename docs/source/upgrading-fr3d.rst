Upgrading FR3D
==============

To properly update the FR3D matlab package, the following steps must be taken.

#. Upgrade the FR3D submodule locally and commit
#. Update the pipeline code on rnatest/rnaprod

Assuming you are in your local copy of the pipeline and on the 'dev' branch,
updating the local copy looks like:

.. code-block:: shell
    $ cd FR3D
    $ git checkout master
    $ git pull
    $ cd ..
    $ git add FR3D
    $ git commit -m 'Updated FR3D to latest'
    $ git push origin dev

From here you then log into remote machine and update it. An example of doing
this on rnatest is. This will update the remote copy.

.. code-block:: shell
    $ ssh rnatest
    $ sudo su - pipeline
    $ cd hub-core
    $ git pull origin dev
    $ git submodule update

It may be useful to examine:

- https://git-scm.com/book/en/v2/Git-Tools-Submodules
- https://git-scm.com/docs/git-submodule
- https://stackoverflow.com/questions/5828324/update-git-submodule-to-latest-commit-on-origin

for information on git submodules and how to work with them.

General Notes
-------------

From here I assume that you are currently logged into the correct machine, your
own, rnatest, or rnaprod, you are the correct user, 'pipeline' in the case of
rnatest and rnaprod and that you are in the pipeline directory. On rnatest and
rnaprod this is `/usr/local/pipeline/hub-core`. An example of ssh'ing in to rnatest
to execute the following commands is:

.. code-block:: shell
    $ ssh username@rnatest
    $ sudo su - pipeline
    password:
    $ cd hub-core

From there all following commands can be run.

Dealing with long running commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some of following commands can take a long time to run. In this case it is
useful to have a `screen` session running if you are on a remote sever. If you
do not have a screen session running then when the ssh connection times out the
program you are running will be killed. If you are not familiar with screen
these links may be helpful.

- https://www.gnu.org/software/screen/manual/screen.html
- http://aperiodic.net/screen/quick_reference
- https://www.linux.com/learn/taking-command-terminal-gnu-screen
- https://en.wikipedia.org/wiki/GNU_Screen
- https://www.linode.com/docs/networking/ssh/using-gnu-screen-to-manage-persistent-terminal-sessions/

Here is how to setup a session rnatest.

.. code-block:: shell
    $ ssh username@rnatest
    $ screen
    $ sudo su - pipeline
    password:
    $ cd hub-core

When you reconnect you can reattach to the screen session using (assuming you
only have 1 session running):

.. code-block:: shell
    $ ssh username@rnatest
    $ screen -r

This will reattach to the session that is logged in as the pipeline user and
you will not have to type in your password. Logging outout to a file
^^^^^^^^^^^^^^^^^^^^^^^^

In several places below I run commands of the form:

.. code-block:: shell
    $ bin/pipeline.py run [options] <command>

These commands will produce *a lot* of output. It is probably useful to log to
a file instead of the screen by doing:

.. code-block:: shell
    $ bin/pipeline.py --log-file some-filename.log --log-mode w run [options] <command>

The ``--log-file some-filename.log`` will log to a file, while the ``--log-mode w``
will overwrite the file when the command is run. By default it appends to an
existing file.

Emailing pipeline logs
^^^^^^^^^^^^^^^^^^^^^^

Also, it is useful to send an email once the command is done. This is done
automatically, but the email will be sent to where it is configured in the
``conf/motifatlas.json`` file. This can be override at the command line using:

.. code-block:: shell
    $ bin/pipeline.py --send-to person@email.com --log-file some-filename.log --log-mode w run [options] <command>

The ``--send-to person@email.com`` will control where the email is sent. It is
a good idea to use this if whomever is running these upgrades is not the person
who normally gets the pipeline emails. The emailing is still a bit experimental
as it can fail when trying to include a log file that is too large for the
setup of BGSU's severs.

Exploring what will be run
^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, there are lots of command that are run below. It is possible, and
sometimes useful, to explore what the pipeline will run with the explore
command. So for example:

.. code-block:: shell
    $ bin/pipeline.py explore  --skip-stage loops.release --recalculate loops.extractor --recalculate loops.positions loops.positions
    download
    pdbs.info
    export.cifatom
    units.info
    mat_files
    loops.extractor  Will Recalculate
    loops.positions  Will Recalculate

The command is a bit slow now, but still useful.

Updating Interactions
---------------------

If the updated FR3D has new interactions then this section has to be done.

There are 3 stages, ``interactions.pairwise``, ``interactions.summary`` and
``interactions.flanking`` which depend on interaction annotations directly.
These all must be recomputed. Below is an example executing each one
independently.

.. code-block:: shell
    $ bin/pipeline.py run --recalculate mat_files --recalculate . --all interactions.pairwise
    $ bin/pipeline.py run --skip-dependencies --recalculate . --all interactions.summary
    $ bin/pipeline.py run --skip-dependencies --recalculate . --all interactions.flanking

Note the ``--recalculate mat_files`` in the first command. This is required
because otherwise the precomputed data that FR3D produces will not be replaced.
If it is not replaced then the same interactions will be imported to the
database for previously competed files. It is very possible to compress the
above steps to a single command like:

.. code-block:: shell
    $ bin/pipeline.py run --recalculate mat_files \
        --recalculate interactions.pairwise \
        --recalculate interactions.summary \
        --recalculate interactions.flanking --all interactions.loader

Once interactions are updated successfully the loops and IFE's will have to be
recomputed.

Updating Loops
--------------

Loops have to be re-extracted if the interaction annotation procedures or the
loop extraction procedures change. This will require updating both the loops
and the loop positions. This can be done in one command like:

.. code-block:: shell
    $ bin/pipeline.py run --skip-stage loops.release --recalculate loops.extractor --recalculate . loops.positions


The command reruns 2 stages, the ``loops.extractor`` which extracts loops from
structures, creates the required mat files and writes to the database, as well
as, the ``loops.positions`` stage which writes the nucleotide to loop
correspondences in the database. This will not remove old loops as we assume
that a loop id is permanent. We skip ``loops.release`` because we are not
trying to create a new loop release, that is done as a regular part of the
update pipeline along with loops quality checks.

Updating IFE's
--------------

IFE's have to be recomputed if the interaction annotations have changed. This
is because IFE's use the number of interactions as part of the process of
building them.

.. code-block:: shell
    $ bin/pipeline.py run --recalculate . ife.info

This part is the trickiest. It is possible that this will through exceptions
for trying to break database constraints. Someone will have to sort those out.
