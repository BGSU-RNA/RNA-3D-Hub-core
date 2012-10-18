"""

A class for downloading .pdb, .cif, and .pdb1 files for all RNA-containing
3D structures.

Example:
python PdbDownloader.py ~/Desktop/f1 ~/Desktop/f2 ~/Desktop/f3

"""

import logging
import os
import sys
import urllib2
import gzip
import shutil
import glob


import PdbInfoLoader
from MotifAtlasBaseClass import MotifAtlasBaseClass


class PdbDownloader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        """
            locations is where pdbs will be placed
            pdbs is an array the files to download
        """
        MotifAtlasBaseClass.__init__(self)
        self.baseurl = 'http://www.rcsb.org/pdb/files/'
        self.filetypes = ['.pdb', '.pdb1']
#         self.filetypes = ['.pdb', '.pdb1', '.cif']
        self.locations = []
        self.pdbs = []
        self.config['email']['subject'] = 'Pdb file sync'

    def get_pdb_list(self):
        """
            Get the latest list of all RNA-containing 3D structures from PDB
        """
        p = PdbInfoLoader.PdbInfoLoader()
        p.get_all_rna_pdbs()
        self.pdbs = p.pdbs
        logging.info('%i RNA 3D structures found in PDB' % len(self.pdbs))

    def set_locations(self, locations):
        """
            Fail if any of the locations doesn't exist
        """
        for location in locations:
            location = os.path.expanduser(location)
            if os.path.exists(location):
                self.locations.append(location)
            else:
                logging.critical("Location %s doesn't exist" % location)
                self.send_report()
                sys.exit(1)
        logging.info('Saving files in %s' % ', '.join(locations))

    def download_files(self):
        """
            Downloads .pdb, .pdb1, and .cif files to the first location, then
            copies all the files over to all the other locations.
        """
        for pdb_id in self.pdbs:
            [self.download(pdb_id, x) for x in self.filetypes]
            if len(self.locations) > 1:
                self.make_copies(pdb_id)

    def download(self, pdb_id, file_type):
        """
            Tries to download the gzipped version of the file. Will crash if
            .pdb or .cif files are not found.
        """
        filename = pdb_id + file_type
        destination = os.path.join(self.locations[0], filename)

        if os.path.exists(destination):
            logging.info('Pdb %s already downloaded' % filename)
            return
        try:
            response = urllib2.urlopen(self.baseurl + filename + '.gz')
        except urllib2.HTTPError:
            if file_type in ['.pdb', '.cif']:
                logging.critical('Pdb file %s could not be downloaded' % pdb_id)
                self.send_report()
                sys.exit(1)
            else:
                logging.info('%s file not found for %s' % (file_type, pdb_id))
                return

        # save gzip
        f = open(destination, 'w')
        f.write(response.read())
        f.close()

        # decompress gzip
        f = gzip.open(destination, 'r')
        decompressed = f.read()
        f.close()

        # save text
        f = open(destination, 'w')
        f.write(decompressed)
        f.close()
        logging.info('Downloaded %s' % destination)

    def make_copies(self, pdb_id):
        """
            Copy files from the first location to all the other locations
        """
        for location in self.locations[1:]:
            source_files = glob.glob(os.path.join(self.locations[0], pdb_id + '*'))
            for source in source_files:
                (head, tail) = os.path.split(source)
                destination = os.path.join(location, tail)
                if os.path.exists(destination):
                    continue
                shutil.copy(source, destination)
            logging.info('Copied %s files to %s' % (pdb_id, location))


def main(argv):
    """
    """
    d = PdbDownloader()
    d.start_logging()
    d.set_locations(argv)
    d.get_pdb_list()
    d.download_files()
    logging.info('Successful update')
    d.send_report()


if __name__ == "__main__":
    main(sys.argv[1:])
