# About
RNA-3D-Hub-core is the backend for the [RNA 3D Hub](http://rna.bgsu.edu/rna3dhub).

RNA-3D-Hub-core contains:

1. [FR3D](https://github.com/BGSU-RNA/FR3D) is included as a submodule.

2. Matlab code for extracting and clustering RNA 3D motifs.

3. Python code responsible for importing non-redundant lists and motif atlas
releases into the database, and id assignment to motifs and non-redundant
equivalence classes.

## Requirements
* python 2.5 or newer (not tested with Python 3)
* matlab R2007b or newer
* MySQL server
* [mlabwap](http://mlabwrap.sourceforge.net) for linking Matlab with Python
* SQLAlchemy for connecting to the database
* [VARNA](http://varna.lri.fr) v3.7 is included to allow for 2D structure generation
* _optional_: nosetests for running Python unit tests

## Installation

1. Download the source code:

    git clone https://github.com/AntonPetrov/RNA-3D-Hub-core.git

2. Create a config file

        pymotifs/motifatlas.cfg

    according to the template found in:

        pymotifs/motifatlas.cfg.txt

3. Create the MySQL database(s) specified in the config file.

4. Submodule initialization

        git submodule init
        git submodule update

## Testing

The software suite includes test datasets for motifs and non-redundant lists.

Unit tests create a special testing environment and shouldn't interfere with
the development or the production versions of the resource. When the unittest
module is imported, the programs will connect to the test database specified
in the config file.

To run unit tests with nosetests:

    nosetests --nologcapture -s pymotifs

* The `--nologcapture` option ensures that all output is logged to a file and is
not intercepted by nosetests. The log file is emailed at the end of the test.
* The `-s` option routes all STDOUT output to the screen for easier monitoring
of test progress.

## Logging and email notifications

Both Python and Matlab programs add their log messages to a file

    MotifAtlas/logs/rna3dhub_log.txt

The log file is refreshed each time the programs are run.

Some programs email this log file using the information specified in the email
section of the config file.

## Directory structure

Write later

## Usage

The main program that triggers the update is `MotifAtlas.py`

When run as a cronjob, need to export a system variable like so:

    export MLABRAW_CMD_STR=/Applications/MATLAB_R2007b/bin/matlab