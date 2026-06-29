# About

RNA-3D-Hub-core is the backend for the [RNA 3D Hub](http://rna.bgsu.edu/rna3dhub).

RNA-3D-Hub-core contains:

1. Code to download, annotate, analyze, and organize RNA- and DNA-containing 3D structures from the PDB each week.  The weekly pipeline run generates equivalence classes of RNA structures representative sets of RNA 3D structures.  As of release 4.0 on 2025-08-13, all code is written in Python 3; for earlier Matlab code see the "master" branch.
2. Every four weeks, code that generates the RNA 3D Motif Atlas is run.  This code was updated considerably shortly before release 4.0.
3. [fr3d-python](https://github.com/BGSU-RNA/fr3d-python) is included as a submodule.  This is an all-python implementation that reads 3D structures, computes symmetry-generated coordinates, and classifies pairwise interactions.

## Documentation

Detailed documentation on this can be found at [readthedocs](http://rna-3d-hub-core.readthedocs.io/).

## Requirements
* python 3.6 or newer
* MySQL server
* SQLAlchemy for connecting to the database
* [fr3d.py](https://github.com/BGSU-RNA/fr3d-python)
* _optional_: py.test for running Python unit tests, but they have not been updated for years

## Installation - very old instructions

1. Download the source code:

    $ git clone https://github.com/AntonPetrov/RNA-3D-Hub-core.git

2. [blank]

3. Create a config file, using the template found in:
   `conf/motifatlas.json.txt`. By default the pipeline will read the config file
in `conf/motifatlas.json`, but this can be changed.

4. Create the MySQL database(s) specified in the config file.

5. Submodule initialization
        git submodule init
        git submodule update

6. Install all dependencies. This can be done by doing
`pip install -f requirements.txt`. This will include nosetests for testing as
well.

## Testing

The software suite includes test datasets for motifs and non-redundant lists.

Unit tests create a special testing environment and shouldn't interfere with the
development or the production versions of the resource. Testing requires that a
`conf/bootstrap.json` config file exists. This will be read for all configuration
information. It is strongly recommended that this be a separate database from
the production database.

To run unit tests with py.test:

```sh
$ py.test
```

## Logging and email notifications

Python logging is configurable. See the section on Configuration for details on
how to define the file to log to. By default logging goes to stdout.

The log file is refreshed each time the programs are run.

## Directory structure

All python code is stored in `pymotifs/`. 
The folder `cif_files` is used to store all downloaded cif files. In
addition, any cached files will be stored there as well. `MotifAtlas` is used to store motif atlas
specific files.

## Usage

The main program that triggers the update is `bin/pipeline`. To run a full
update, do `bin/pipeline update --all`. This will run all stages of the pipeline
on all RNA containing structures. If you want to only run one stage of the
pipeline, such as downloading all pdbs, simply do `bin/pipeline download --all`.
For details on what stages are available do `bin/pipeline help`.
If you only want to run the pipeline on one file or more files, simply do:
`bin/pipeline download 2AW7 1J5E`.

When run as a cronjob, you must export a system variable that defines the matlab
command string for mlabwrap like so:

    export MLABRAW_CMD_STR=/Applications/MATLAB_R2007b/bin/matlab

An example crontab file that runs the update every friday can be found in
`files/examples/update.cron`.

## Configuration

This requires that any firewall you are using allow passive FTP connections.

## License

Source code: MIT License
Documentation and published data: Creative Commons Attribution 4.0 International (CC BY 4.0)
