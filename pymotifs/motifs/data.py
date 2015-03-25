"""

Main script for loading motif groups into RNA 3D Hub.
Manages id tracking, version numbers etc.
Treats internal, hairpin and junction loops separately.

The folders must be prefixed with the loop type,
for example: il1, il2, hl1, hl2. The release precedence is determined by
sorting folder names.

Usage: python MotifLoader [options]

Options:
  -t                    motif type (mandatory): il, hl, j3, j4 etc
  -h, --help            show this help

Examples:
python MotifLoader.py -t il
python MotifLoader.py -t hl

"""

import os
import fnmatch
import logging

from sqlalchemy import desc

from CollectionsMerger import CollectionsMerger as MotifCollectionMerger
from pymotifs.motifs.collections import MotifCollection
from pymotifs.motifs.uploader import Uploader

from pymotifs.models import MlReleases

logger = logging.getLogger(__name__)


class MotifLoader(object):

    def __init__(self, motif_type=''):
        self.success = False
        self.motifs_root = self.config['locations']['releases_dir']
        self.motif_type = motif_type.upper()
        self.done = []
        self.folders = []

    def __get_data_folders(self):
        """folders must begin with loop_type (IL or HL)"""
        for file in os.listdir(self.motifs_root):
            if fnmatch.fnmatch(file, self.motif_type + '*'):
                self.folders.append(file)
                logger.info(file)
        self.folders = sorted(self.folders)

    def __initialize_files(self, folder):
        """
            This function generates full paths to all files that make up a
            motif atlas release. These files are produced by Matlab are
            are expected to be found in each release directory.
        """
        self.f = dict()
        self.f['description'] = folder
        filenames = {
            'MotifList':         'MotifList.csv',
            'MotifLoopOrder':    'MotifLoopOrder.csv',
            'MotifPositions':    'MotifPositions.csv',
            'MutualDiscrepancy': 'MutualDiscrepancy.csv',
            'BpSignatures':      'MotifBpSignatures.csv',
            'MatFiles_origin':   'Groups',
            'MatFiles_dest':     'mat',
            '2ds_origin':        '2ds',
            'correspondences':   'correspondences.txt'
        }

        self.f['folder'] = os.path.join(self.motifs_root, folder)
        for k, v in filenames.iteritems():
            self.f[k] = os.path.join(self.f['folder'], v)

        # location for images accessible to the webserver
        self.f['2ds_destination'] = self.config['locations']['2ds_destination']
        if not os.path.exists(self.f['2ds_destination']):
            os.mkdir(self.f['2ds_destination'])

    def get_processed_releases(self):
        """The `description` field corresponds to the folder name that the data
        for the release came from
        """

        with self.session() as session:
            query = session.query(MlReleases).\
                filter(MlReleases.type == self.motif_type)

        for release in query:
            self.done.append(release.description)

    def list_all_releases(self, type):
        all_releases = []
        with self.session() as session:
            query = session.query(MlReleases).\
                filter(MlReleases.type == type).\
                order_by(desc(MlReleases.date))

            for release in query:
                all_releases.append(release.id)
        return all_releases

    def import_motif_release(self):

        c1 = MotifCollection(file=self.f['MotifList'], type=self.motif_type)
        c2 = MotifCollection(release='latest', type=self.motif_type)

        uploader = Uploader(ensembles=MotifCollectionMerger(c1, c2),
                            mode=self.config['release_mode']['motifs'],
                            description=self.f['description'],
                            files=self.f,
                            motif_type=self.motif_type)

        return uploader.data()

    # def compare_all_releases(self):
    #     all_releases = self.list_all_releases(self.motif_type)
    #     if len(all_releases) <= 2:
    #         return

    #     all_releases = all_releases[2:]
    #     c1 = MotifCollection(release='latest', type=self.motif_type)
    #     for release in all_releases:
    #         c2 = MotifCollection(release=release, type=self.motif_type)
    #         logger.info('Comparing %s and %s' % (c1.release, c2.release))
    #         A = Uploader(ensembles=MotifCollectionMerger(c1, c2),
    #                      upload_mode='release_diff',
    #                      motif_type=self.motif_type)
    #         A.import_release()

    def import_data(self):
        self.get_processed_releases()
        self.__get_data_folders()
        for folder in self.folders:
            if len(folder) < 6:
                logger.info('Skipping folder %s' % folder)
                continue

            if folder in self.done:
                logger.info('Already imported %s' % folder)
                continue

            self.__initialize_files(folder)
            self.import_motif_release()
            # self.compare_all_releases()
        self.success = True


def main(argv):

    opts = []
    motif_type = None
    for opt, arg in opts:
        if opt in ('-t', '--type'):
            motif_type = arg

    L = MotifLoader(motif_type=motif_type)
    L.import_data()
