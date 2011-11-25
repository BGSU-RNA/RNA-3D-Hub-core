"""Motifs loader
Loads motif groups into a database. Manages id tracking, version numbers etc.

Usage: python aMLMotifLoader [options]

Options:
  -d                    drop all tables and reload
  -h, --help            show this help

Examples:
python aMLMotifLoader.py
"""

import pdb, sys, getopt, os, re
from aMLSqlAlchemyClasses import *

__author__ = 'Anton Petrov'



class Loader:
    """
    """

    def __init__(self):
        self.motifs_root = '/FR3D/Workspace/Releases'

#         self.done  = list_done()
        self.folders = sorted(os.listdir(self.motifs_root))


    def initialize_files(self, folder):
        """
        """
        self.f = dict()
        self.f['folder']            = os.path.join(self.motifs_root, folder) #'/FR3D/Workspace/Releases/iljun6'
        self.f['MotifList']         = os.path.join(self.f['folder'],'MotifList.csv')
        self.f['MotifLoopOrder']    = os.path.join(self.f['folder'],'MotifLoopOrder.csv')
        self.f['MotifPositions']    = os.path.join(self.f['folder'],'MotifPositions.csv')
        self.f['MutualDiscrepancy'] = os.path.join(self.f['folder'],'MutualDiscrepancy.csv')
        self.f['MatFiles']          = os.path.join(self.f['folder'],'mat')


    def import_motif_release(self):
        """
        """
#                     A = Uploader(ensembles=NRCollectionMerger(c1,c2),
#                                  release_mode='reuse',
#                                  release_description=folder)
#       c1 = NR_eqclass_collection(file=self.temp_file, resolution=self.resolution_labels[res_id])
#       c2 = NR_eqclass_collection(release='previous',resolution=self.resolution_labels[res_id])

        c1 = MotifCollection(file=self.f['MotifList'])
        c2 = MotifCollection(release='latest')
        pdb.set_trace()
        A = Uploader(ensembles=MotifCollectionMerger(c1,c2),
                        mode='minor',
                        description=self.f['folder'],
                        files=self.f)


    def import_data(self):
        """
        """
        for folder in self.folders:
        	if len(folder) == 6:
				self.initialize_files(folder)
				self.import_motif_release()



def usage():
    print __doc__


def main(argv):
    """
    """
    try:
        opts, args = getopt.getopt(argv, "d", ['help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-d':
            print 'Confirm dropping all tables and reloading the data by pressing c'
            pdb.set_trace()
            drop_all()
            Base.metadata.create_all(engine)
        elif opt in ('-h', '--help'):
            usage()
            sys.exit()


    L = Loader()
    L.import_data()


if __name__ == "__main__":
    main(sys.argv[1:])

