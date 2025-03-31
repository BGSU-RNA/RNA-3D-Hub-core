"""
Download CIF files.

This will download compressed cif files and place under PDBFiles in the defined
FR3D directory.
"""

from contextlib import contextmanager
import os
import urllib.request

from pymotifs import core
from pymotifs import utils


class Writer(core.FileHandleSaver):
    @contextmanager
    def writer(self, pdb, **kwargs):
        with super(Writer, self).writer(pdb, **kwargs) as raw:
            yield raw.write


class Downloader(core.Loader):
    name = 'downloader'
    file_url = 'https://files.rcsb.org/download/{pdb}.cif.gz'
    update_gap = False
    dependencies = set()
    saver = Writer

    def __init__(self, *args, **kwargs):
        super(Downloader, self).__init__(*args, **kwargs)
        self.FH = utils.FetchHelper(allow_fail=True)
        self.location = self.config['locations']['cif_files']

    def filename(self, name, **kwargs):
        return os.path.realpath(os.path.normpath(os.path.join(self.location, name + '.cif.gz')))

    def url(self, name, **kwargs):
        return self.file_url.format(pdb=name)

    def remove(self, entry, **kwargs):
        if self.has_data(entry) and not kwargs.get('dry_run'):
            os.remove(self.filename(entry))

    def has_data(self, entry, **kwargs):
        return os.path.exists(self.filename(entry))

    def data(self, name, **kwargs):
        downloaded = False
        try:
            # skip the helpers and saver stage and just get it done
            url = self.url(name, **kwargs)
            filename = self.filename(name, **kwargs)
            self.logger.info('Downloading %s from %s' % (filename,url))
            # next line failed on 2/18/2025, timeout
            urllib.request.urlretrieve(url, filename)
            downloaded = True
        except Exception as err:
            self.logger.error('%s could not be downloaded', name)
            self.logger.exception(err)
            raise core.Skip("Couldn't get %s" % name)

        if downloaded:
            raise core.Skip("Downloaded to %s" % filename)

        # if not content:
        #     core.Skip("Downloaded empty file %s" % name)

        # return content
