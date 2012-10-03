<h3>About</h3>
PyFR3DMotifs is the backend for the <a href="http://rna.bgsu.edu/rna3dhub">RNA 3D Hub</a>.

It contains:

1. core FR3D modules written in Matlab. In the near future PyFR3DMotifs will be integrated with the main <a href="https://github.com/BGSU-RNA/FR3D">FR3D Github repo</a>.

2. additional Matlab code for extraction and clustering of RNA 3D motifs.

3. Python code responsible for importing non-redundant lists and motif atlas releases into the database, and id assignment to motifs and non-redundant equivalence classes.

<h3>Requirements</h3>
* python 2.5 or newer (not tested with Python 3)
* matlab R2007b or newer
* MySQL server
* <a href="http://mlabwrap.sourceforge.net/">mlabwap</a> for linking Matlab with Python
* SQLAlchemy for connecting to the database
* <a href="http://varna.lri.fr">VARNA</a> v3.7 is included to allow for 2D structure generation
* <i>optional</i>: nosetests for running Python unit tests

<h3>Installation</h3>

1. Download the source code:

    git clone https://github.com/AntonPetrov/PyFR3DMotifs.git

1. Create a config file

    pymotifs/motifatlas.cfg

according to the template found in:

    pymotifs/motifatlas.cfg.txt

2. Run unit tests:

cd pymotifs

    nosetests --nologcapture

The software suite includes test datasets for motifs and non-redundant lists.

Unit tests create a special testing environment and shouldn't interfere with the development or the production versions of the resource. When the unittest module is imported, the programs will connect to the test database specified in the config file.

<h3>Usage</h3>

The main program that triggers the update is MotifAtlas.py