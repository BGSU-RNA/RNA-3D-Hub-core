"""
Motifs loader
Loads motif groups into a database. Manages id tracking, version numbers etc.
Treats internal, hairpin and junction loops separately. The root folder must
contain folders with mat files. The folders must be prefixed with the loop type,
for example: il1, il2, hl1, hl2. The release precedence is determined by sorting
folder names.

Usage: python aMLMotifLoader [options]

Options:
  -t                    motif type (mandatory): il, hl, j3, j4 etc
  -d                    drop all tables and reload
  -h, --help            show this help

Examples:
python aMLMotifLoader.py
"""

import pdb, sys, getopt, os, re, fnmatch
from aMLSqlAlchemyClasses import *

__author__ = 'Anton Petrov'



class Loader:
    """
    """

    def __init__(self, motif_type=''):
        self.motifs_root = '/FR3D/Workspace/Releases'
        self.motif_type  = motif_type

        self.done  = list_done(self.motif_type)

        self.folders = []
        for file in os.listdir(self.motifs_root):
            if fnmatch.fnmatch(file, self.motif_type + '*'):
                self.folders.append(file)
                print file
        self.folders = sorted(self.folders)


    def initialize_files(self, folder):
        """
        """
        self.f = dict()
        self.f['description']       = folder
        self.f['folder']            = os.path.join(self.motifs_root, folder) #'/FR3D/Workspace/Releases/iljun6'
        self.f['MotifList']         = os.path.join(self.f['folder'],'MotifList.csv')
        self.f['MotifLoopOrder']    = os.path.join(self.f['folder'],'MotifLoopOrder.csv')
        self.f['MotifPositions']    = os.path.join(self.f['folder'],'MotifPositions.csv')
        self.f['MutualDiscrepancy'] = os.path.join(self.f['folder'],'MutualDiscrepancy.csv')
        self.f['MatFiles']          = os.path.join(self.f['folder'],'mat')
        self.f['2ds_origin']        = os.path.join(self.f['folder'],'html','2ds')
        self.f['2ds_destination']   = '/Servers/rna.bgsu.edu/img/MotifAtlas'


    def import_motif_release(self):
        """
        """
        c1 = MotifCollection(file=self.f['MotifList'],type=self.motif_type)
        c2 = MotifCollection(release='latest',type=self.motif_type)
        A = Uploader(ensembles=MotifCollectionMerger(c1,c2),
                        mode='minor',
                        description=self.f['description'],
                        files=self.f,
                        motif_type=self.motif_type)


    def compare_all_releases(self):
        """
        """
        all_releases = list_all_releases(self.motif_type)
        if len(all_releases) <= 2:
            return
        all_releases = all_releases[2:]
        c1 = MotifCollection(release='latest',type=self.motif_type)
        for release in all_releases:
            c2 = MotifCollection(release=release,type=self.motif_type)
            print c1.release, c2.release
            A = Uploader(ensembles=MotifCollectionMerger(c1,c2),
                         upload_mode='release_diff',
                         motif_type=self.motif_type)


    def import_data(self):
        """
        """
        for folder in self.folders:
            if len(folder) != 6:
                print 'Skipping folder ', folder
                continue

            if folder in self.done:
                print 'Already imported ', folder
                continue

            self.initialize_files(folder)
            self.import_motif_release()
            self.compare_all_releases()



def usage():
    print __doc__


def main(argv):
    """
    """
    try:
        opts, args = getopt.getopt(argv, "dht:", ['help','type'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    motif_type = None
    for opt, arg in opts:
        if opt == '-d':
            print 'Confirm dropping all tables and reloading the data by pressing c'
            pdb.set_trace()
            drop_all()
            Base.metadata.create_all(engine)
        elif opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-t', '--type'):
            motif_type = arg

    if motif_type is None:
        usage()
        print 'Please specify motif type: il, hl, j3, j4 etc'
        sys.exit(2)


    L = Loader(motif_type=motif_type)
    L.import_data()
#     drop_temp_tables()
#     pdb.set_trace()

if __name__ == "__main__":
    main(sys.argv[1:])

