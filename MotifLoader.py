"""
Motif loader

Loads motif groups into the database. Manages id tracking, version numbers etc.
Treats internal, hairpin and junction loops separately. The root folder must
contain folders with mat files.

The folders must be prefixed with the loop type,
for example: il1, il2, hl1, hl2. The release precedence is determined by sorting
folder names.

Usage: python MotifLoader [options]

Options:
  -t                    motif type (mandatory): il, hl, j3, j4 etc
  -h, --help            show this help

Examples:
python MotifLoader.py
"""

import pdb
import sys
import getopt
import os
import re
import fnmatch
import logging

from sqlalchemy import desc


from MLSqlAlchemyClasses import session, MotifCollection, Uploader, Release, \
                                MotifCollectionMerger
from MotifAtlasBaseClass import MotifAtlasBaseClass


__author__ = 'Anton Petrov'


class MotifLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self, motif_type=''):
        MotifAtlasBaseClass.__init__(self)
        self.motifs_root = self.config['locations']['releases_dir']
        self.motif_type  = motif_type.upper()
        self.done  = self.get_processed_releases(self.motif_type)
        self.folders = []
        self.__get_data_folders()

    def __get_data_folders(self):
        """folders must begin with loop_type (IL or HL)"""
        for file in os.listdir(self.motifs_root):
            if fnmatch.fnmatch(file, self.motif_type + '*'):
                self.folders.append(file)
                print file
        self.folders = sorted(self.folders)

    def __initialize_files(self, folder):
        """
        """
        self.f = dict()
        self.f['description']       = folder
        self.f['folder']            = os.path.join(self.motifs_root, folder)
        self.f['MotifList']         = os.path.join(self.f['folder'],'MotifList.csv')
        self.f['MotifLoopOrder']    = os.path.join(self.f['folder'],'MotifLoopOrder.csv')
        self.f['MotifPositions']    = os.path.join(self.f['folder'],'MotifPositions.csv')
        self.f['MutualDiscrepancy'] = os.path.join(self.f['folder'],'MutualDiscrepancy.csv')
        self.f['BpSignatures']      = os.path.join(self.f['folder'],'MotifBpSignatures.csv')
        self.f['MatFiles_origin']   = os.path.join(self.f['folder'],'Groups')
        self.f['MatFiles_dest']     = os.path.join(self.f['folder'],'mat')
        self.f['2ds_origin']        = os.path.join(self.f['folder'],'2ds')
        self.f['2ds_destination']   = '/Servers/rna.bgsu.edu/img/MotifAtlas'
        self.f['correspondences']   = os.path.join(self.f['folder'],'correspondences.txt')

    def get_processed_releases(self, type):
        """the `description` field corresponds to the folder name that the data
        for the release came from"""
        done = []
        for release in session.query(Release).filter(Release.type==type).all():
            done.append(release.description)
        return done

    def list_all_releases(self, type):
        all_releases = []
        for release in session.query(Release). \
                               filter(Release.type==type). \
                               order_by(desc(Release.date)). \
                               all():
            all_releases.append(release.id)
        return all_releases

    def import_motif_release(self):
        """
        """
        c1 = MotifCollection(file=self.f['MotifList'], type=self.motif_type)
        c2 = MotifCollection(release='latest', type=self.motif_type)
        A = Uploader(ensembles=MotifCollectionMerger(c1,c2),
                        mode=self.config['release_mode']['motifs'], #'minor',
                        description=self.f['description'],
                        files=self.f,
                        motif_type=self.motif_type)
        A.import_release()

    def compare_all_releases(self):
        """
        """
        all_releases = self.list_all_releases(self.motif_type)
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
            A.import_release()

    def import_data(self):
        """
        """
        for folder in self.folders:
            if len(folder) < 6:
                print 'Skipping folder ', folder
                continue

            if folder in self.done:
                print 'Already imported ', folder
                continue

            self.__initialize_files(folder)
            self.import_motif_release()
            self.compare_all_releases()


def usage():
    print __doc__


def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    try:
        opts, args = getopt.getopt(argv, "t:r:h", ['help','type'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    motif_type = None
    release_to_remove = None
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-t', '--type'):
            motif_type = arg
        elif opt in ('-r'):
            release_to_remove = arg

    if motif_type is None:
        usage()
        print 'Please specify motif type: il, hl, j3, j4 etc'
        sys.exit(2)

    if release_to_remove is not None:
        print 'Are you sure you want to delete %s release %s?. Press c to continue or q to abort.' \
              % (motif_type, release_to_remove)
        pdb.set_trace()
        U = Uploader(motif_type=motif_type)
        U.remove_release(release_to_remove)
        print 'Success'
        return

    L = MotifLoader(motif_type=motif_type)
    L.import_data()


if __name__ == "__main__":
    main(sys.argv[1:])

