"""
A program for downloading .cif, files for all RNA-containing
3D structures.
"""

import os

import core
import utils


class Downloader(core.Loader):
    file_url = 'http://www.rcsb.org/pdb/files/'
    name = 'downloader'
    update_gap = False

    def __init__(self, config, maker):
        self.gzip = utils.GzipFetchHelper(allow_fail=True)
        self.location = os.path.join(config['locations']['fr3d_root'],
                                     'PDBFiles')
        super(core.Loader, self).__init__(config, maker)

    def filename(self, name):
        return os.path.join(self.location, name + '.cif')

    def url(self, name):
        return self.file_url + name + '.cif.gz'

    def remove(self, entry):
        if self.has_data(entry):
            os.remove(self.filename(entry))

    def has_data(self, entry):
        return os.path.exists(self.filename(entry))

    def store(self, data):
        filename, text = data
        with open(filename, 'w') as out:
            out.write(text)
        self.logger.info('Downloaded %s' % filename)

    def data(self, name, **kwargs):
        destination = self.filename(name)

        try:
            content = self.gzip(self.url(name))
        except:
            self.logger.error('%s could not be downloaded', name)
            raise core.SkipValue("Couldn't get %s" % name)

        if not content:
            core.SkipValue("Downloaded empty file %s" % name)

        return destination, content
