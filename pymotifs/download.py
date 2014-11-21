"""
A program for downloading .pdb, .cif, and .pdb1 files for all RNA-containing
3D structures.
"""

import os
import xml.dom.minidom

import core
import utils


class Downloader(core.Loader):
    file_url = 'http://www.rcsb.org/pdb/files/'
    entity_info_url = 'http://www.pdb.org/pdb/rest/getEntityInfo?structureId='
    name = 'downloader'
    update_gap = False

    def __init__(self, config, maker):
        self.helper = utils.WebRequestHelper(parser=self.parse)
        self.gzip = utils.GzipFetchHelper(allow_fail=True)
        self.location = os.path.join(config['locations']['fr3d_root'],
                                     'PDBFiles')
        super(core.Loader, self).__init__(config, maker)

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
            self.logger.error('Bioassembly query failed for  %s' % pdb)
            return None

    def remove(self, entry):
        pass

    def transform(self, pdb):
        values = [(pdb, '.cif'), (pdb, '.pdb')]
        if self.assembly_count(pdb):
            values.append((pdb, '.pdb1'))
        return values

    def has_data(self, filename):
        full = os.path.join(self.location, filename)
        return os.path.exists(full)

    def store(self, data):
        filename, text = data
        with open(filename, 'w') as out:
            out.write(text)
        self.logger.info('Downloaded %s' % filename)

    def data(self, filename, **kwargs):
        pdb, extension = filename
        destination = os.path.join(self.location, pdb + extension)

        try:
            content = self.gzip(self.file_url + pdb + extension + '.gz')
        except:
            if extension == '.pdb' or extension == '.cif':
                self.logger.critical('%s could not be downloaded',
                                     pdb + extension)
                raise

        if not content:
            core.SkipValue("Downloaded empty file %s" % pdb + extension)

        return destination, content
