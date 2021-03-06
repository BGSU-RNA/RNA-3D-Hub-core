"""Download CIF files.

This will download compressed cif files and place under PDBFiles in the defined
FR3D directory.
"""

import os
from contextlib import contextmanager

from pymotifs import core
from pymotifs import utils


class Writer(core.FileHandleSaver):
    @contextmanager
    def writer(self, pdb, **kwargs):
        with super(Writer, self).writer(pdb, **kwargs) as raw:
            yield raw.write


class Downloader(core.Loader):
    file_url = 'http://files.rcsb.org/download/{pdb}.cif.gz'
    name = 'downloader'
    update_gap = False
    dependencies = set()
    saver = Writer

    def __init__(self, *args, **kwargs):
        super(Downloader, self).__init__(*args, **kwargs)
        self.gzip = utils.GzipFetchHelper(allow_fail=True)
        self.location = os.path.join(self.config['locations']['fr3d_root'],
                                     'PDBFiles')

    def filename(self, name, **kwargs):
        return os.path.realpath(os.path.normpath(os.path.join(self.location, name + '.cif')))

    def url(self, name, **kwargs):
        return self.file_url.format(pdb=name)

    def remove(self, entry, **kwargs):
        if self.has_data(entry) and not kwargs.get('dry_run'):
            os.remove(self.filename(entry))

    def has_data(self, entry, **kwargs):
        return os.path.exists(self.filename(entry))

    def data(self, name, **kwargs):
        try:
            content = self.gzip(self.url(name, **kwargs))
        except Exception as err:
            self.logger.error('%s could not be downloaded', name)
            self.logger.exception(err)
            raise core.Skip("Couldn't get %s" % name)

        if not content:
            core.Skip("Downloaded empty file %s" % name)

        return content
