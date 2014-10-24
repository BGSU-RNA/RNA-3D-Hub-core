"""

A program for downloading .pdb, .cif, and .pdb1 files for all RNA-containing
3D structures.

Downloads the files specified in self.filetypes to the specified directory.

Example:
python PdbDownloader.py ~/Desktop/f1

"""

import os
import logging
import xml.dom.minidom

import utils

logger = logging.getLogger(__name__)


class Downloader(object):
    file_url = 'http://www.rcsb.org/pdb/files/'
    entity_info_url = 'http://www.pdb.org/pdb/rest/getEntityInfo?structureId='

    def __init__(self, config, session):
        """
            locations is where pdbs will be placed
            pdbs is an array the files to download
        """
        self.config = config
        self.filetypes = ['.pdb', '.pdb1', '.cif']
        self.helper = utils.WebRequestHelper(parser=self.parse)
        self.location = os.path.join(config['locations']['fr3d_root'],
                                     'PDBFiles')

    def parse(self, response):
        dom = xml.dom.minidom.parseString(response.text)
        tag = dom.getElementsByTagName('PDB')[0]
        raw = tag.attributes['bioAssemblies'].value
        return int(raw)

    def assemblies_count(self, pdb):
        """Find the number of biological assemblies associated with a pdb id.
        """
        try:
            return self.helper(self.entity_info_url + pdb)
        except:
            logger.error('Bioassembly query failed for  %s' % pdb)
            return None

    def download(self, pdb, file_type):
        """
            Tries to download the gzipped version of the file. Will crash if
            .pdb or .cif files are not found.
        """
        filename = pdb + file_type

        if file_type == '.pdb1':
            count = self.assemblies_count(pdb)
            if count == 0:
                return

            if count is None:
                logger.warning("Skipping attempt to get %s", filename)
                return

        destination = os.path.join(self.location, filename)

        if os.path.exists(destination):
            logger.info('%s already downloaded' % filename)
            return

        try:
            helper = utils.GzipFetchHelper(allow_fail=True)
            content = helper(self.file_url + filename + '.gz')
            if not content:
                logger.error("Downloaded an empty file for %s", filename)
                return

            with open(destination, 'w') as out:
                out.write(content)
        except:
            if file_type == '.pdb' or file_type == '.cif':
                logger.critical('%s could not be downloaded', filename)
                raise
            else:
                logger.warning('%s not found', filename)
                return

        logger.info('Downloaded %s' % destination)

    def __call__(self, pdbs, **kwargs):
        """
            Downloads .pdb, .pdb1, and .cif files to the first location, then
            copies all the files over to all the other locations.
        """
        if not os.path.isdir(self.location):
            logger.critical("%s does not exist", self.location)
            raise Exception("Can't download to missing directory")

        for pdb in pdbs:
            for filetype in self.filetypes:
                try:
                    self.download(pdb, filetype)
                except:
                    logger.error("Failed to download %s", pdb)
                    continue
