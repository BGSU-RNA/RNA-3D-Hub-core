"""

NR list loader
Loads NR groups into a database. Manages id tracking, version numbers etc.

Usage: python NRClassLoader [options]

Options:
  -d                    drop all tables and reload
  -h, --help            show this help

Examples:
python NRClassLoader.py

"""

import pdb
import sys
import getopt
import os
import re
import logging

from nratlas.NRSqlAlchemyClasses import *
from nratlas.NRUploader import Uploader
from nratlas.NRCollections import NR_eqclass_collection, NRCollectionMerger
from MotifAtlasBaseClass import MotifAtlasBaseClass

__author__ = 'Anton Petrov'


class Loader(MotifAtlasBaseClass):
    """
    """

    def __init__(self):
        """
        """
        MotifAtlasBaseClass.__init__(self)
        self.temp_file          = 'temp.csv'
        self.nrlists_root       = self.config['locations']['nrlists_dir']
        self.resolutions        = ['1,5A','2A','2,5A','3A','3,5A','4A','20A','All_Resolution']
        self.resolution_labels  = ['1.5','2.0','2.5','3.0','3.5','4.0','20.0','all']
        self.done = []
        self.list_done()
        self.lists = sorted(os.listdir(self.nrlists_root))
        self.success = False # status of the current update

    def list_done(self):
        """
        """
        for release in session.query(NR_release).all():
            self.done.append(release.description)

    def list_all_releases(self):
        """
        """
        all_releases = []
        for release in session.query(NR_release).\
                               order_by(desc(NR_release.date)).\
                               all():
            all_releases.append(release.id)
        return all_releases

    def parse_file(self, file=''):
        """
        """
        f = open(file, 'r')
        ofn = open(self.temp_file, 'w')
        for (i,line) in enumerate(f):
            tds = []
            if line.find('<tr') != -1:
                tds = line.split('td')
                if len(tds)>1:
                    # parse first td
                    ind = tds[1].find('pdbId=')
                    rep = tds[1][ind+6:ind+10]
                    ofn.write( '%s,Group_%i\n' % (rep, i) )
                    # parse the last td
                    s = tds[-2]
                    if s.find('&nbsp;') != -1:
                        pass
                    elif s.find('<a') == -1:
                        # older format
                        pdbs = re.findall('\w{4}',s)
                        for p in pdbs:
                            ofn.write( '%s,Group_%i\n' % (p, i) )
                    else:
                        # newer table format with links
                        anchors = s.split('<a')
                        for link in anchors:
                            p = re.findall('>(\w{4})<',link)
                            if len(p) > 0 and len(p[0]) == 4: # "Compare" link
                                ofn.write( '%s,Group_%i\n' % (p[0], i) )

    def compare_all_releases(self):
        """
        """
        all_releases = self.list_all_releases()
        if len(all_releases) <= 2:
            return
        all_releases = all_releases[2:]
        for resolution in self.resolution_labels:
            for release in all_releases:
                c1 = NR_eqclass_collection(release='latest', resolution=resolution)
                c2 = NR_eqclass_collection(release=release,  resolution=resolution)
                logging.info('%s, %s, %s' % (c1.release, c2.release, resolution))
                A = Uploader(ensembles=NRCollectionMerger(c1,c2), upload_mode='release_diff')
                A.update_database()

    def import_lists(self,folder):
        """
        """
        for (res_id,resolution) in enumerate(self.resolutions):
            nr_html = 'Nonredundant_' + resolution + '.html'
            full_nr_html = os.path.join(self.nrlists_root, folder, nr_html)
            if os.path.isfile(full_nr_html):
                print full_nr_html
                self.parse_file(file=full_nr_html) # creates self.temp_file

                c1 = NR_eqclass_collection(file=self.temp_file,
                                           resolution=self.resolution_labels[res_id])
                if res_id == 0:
                    c2 = NR_eqclass_collection(release='latest',resolution=self.resolution_labels[res_id])
                    A = Uploader(ensembles=NRCollectionMerger(c1,c2),
                                 release_mode='minor',
                                 release_description=folder)
                    A.update_database()
                else:
                    c2 = NR_eqclass_collection(release='previous',resolution=self.resolution_labels[res_id])
                    A = Uploader(ensembles=NRCollectionMerger(c1,c2),
                                 release_mode='reuse',
                                 release_description=folder)
                    A.update_database()

    def import_data(self):
        """
        """
        for folder in self.lists:
            if len(folder) != 8:
                print 'Skipping folder ', folder
                continue

            if folder in self.done:
                print 'Already imported ', folder
                continue

            self.import_lists(folder)
            self.compare_all_releases()

        if os.path.exists(self.temp_file):
            os.remove(self.temp_file)

        self.success = True
        print 'Successful update'


def usage():
    print __doc__


def main(argv):
    """
    """

    try:
        opts, args = getopt.getopt(argv, "dcr:", ['help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-d':
            print 'Confirm dropping all tables and reloading the data by pressing c or press q to abort:'
            pdb.set_trace()
            drop_all()
            Base.metadata.create_all(engine)
        elif opt == '-r':
            print 'Confirm removing NR release %s by pressing c or press q to abort:' % arg
            pdb.set_trace()
            U = Uploader()
            U.start_logging()
            U.remove_release(arg)
            logging.info('Release removed')
            sys.exit()
        else:
            usage()
            sys.exit()

    L = Loader()
    L.start_logging()
    L.import_data()
    logging.info( '%s completed' % __file__ )
    L.send_report()


if __name__ == "__main__":
    main(sys.argv[1:])
