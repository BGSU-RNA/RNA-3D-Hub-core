"""
A program for downloading .cif files for all RNA-containing
3D structures.
"""

import os

from pymotifs import core
from pymotifs import utils


class Downloader(core.Loader):
    file_url = 'http://www.rcsb.org/pdb/files/'
    name = 'downloader'
    update_gap = False
    dependencies = set()

    def __init__(self, *args, **kwargs):
        super(Downloader, self).__init__(*args, **kwargs)
        self.gzip = utils.GzipFetchHelper(allow_fail=True)
        self.location = os.path.join(self.config['locations']['fr3d_root'],
                                     'PDBFiles')

    def filename(self, name):
        return os.path.join(self.location, name + '.cif')

    def url(self, name):
        return self.file_url + name + '.cif.gz'

    def remove(self, entry, **kwargs):
        if self.has_data(entry):
            os.remove(self.filename(entry))

    def has_data(self, entry, **kwargs):
        return os.path.exists(self.filename(entry))

    def store(self, data, **kwargs):
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
            raise core.Skip("Couldn't get %s" % name)

        if not content:
            core.Skip("Downloaded empty file %s" % name)

        return destination, content
