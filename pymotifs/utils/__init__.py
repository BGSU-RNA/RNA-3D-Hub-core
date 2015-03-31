import os
import logging
from ftplib import FTP
import cStringIO as sio
import itertools as it
import collections as coll
import gzip

import requests

logger = logging.getLogger(__name__)


class MissingFileException(Exception):
    """This is raised when we can't find a file. For example a cif file for a
    pdb does not exist.
    """
    pass


class WebRequestFailed(Exception):
    """Raised when we have failed all attempts at getting a url.
    """


class EmptyResponse(Exception):
    """Raised when processing an empty response.
    """


class RetryFailedException(Exception):
    """Raised when all attempts at retrying something have failed.
    """


class GetAllRnaPdbsError(Exception):
    """Raised with we cannot get all RNA containing PDBS
    """
    pass


class GetCustomReportError(Exception):
    """Raised when we cannot get a custom report
    """
    pass


def grouper(n, iterable):
    """Group an iterable in chunks of a max size.

    :n: The max size.
    :iterable: The iterable to group.
    :returns: A iterable with at most n entries per chunk.
    """

    iterator = iter(iterable)
    while True:
        chunk = tuple(it.islice(iterator, n))
        if not chunk:
            return
        yield chunk


def row2dict(row):
    d = {}
    for column in row.__table__.columns:
        d[column.name] = str(getattr(row, column.name))
    return d


def known(config, pdb=True, cif=True, pdb1=False):
    path = os.path.join(config['locations']['fr3d_root'], 'FR3D', 'PDBFiles')
    names = coll.defaultdict(dict)

    for filename in os.listdir(path):
        if not os.path.isfile(os.path.join(path, filename)):
            continue
        name, ext = os.path.splitext(filename)
        ext = ext.replace('.', '')
        names[name][ext] = True

    for name, exts in names.items():
        if not pdb or (pdb and exts.get('pdb')):
            if not cif or (cif and exts.get('cif')):
                if not pdb1 or (pdb1 and exts.get('pdb1')):
                    yield name


class RetryHelper(object):
    """A base class for retrying an action.
    """

    def __init__(self, retries=3, allow_fail=False, log_last_only=True):
        self.retries = retries
        self.allow_fail = allow_fail
        self.log_last_only = log_last_only

    def __call__(self, *args, **kwargs):
        for index in xrange(self.retries):
            try:
                return self.action(*args, **kwargs)
            except Exception as err:
                if not self.log_last_only or \
                        (self.log_last_only and index == (self.retries - 1)):

                    if not self.allow_fail:
                        logger.warning("Failed retry attempt #%s",
                                       str(index + 1))
                        logger.exception(err)

        if not self.allow_fail:
            logger.error("All attempts at retrying failed")
            raise RetryFailedException()


class WebRequestHelper(RetryHelper):
    """A class to help with making web requests. This deals with the retrying
    and making sure the response body is not empty. If the max number of
    retries is reached we raise an exception. In addition, all steps are
    logged.
    """

    def __init__(self, allow_empty=False, method='get', parser=None, **kwargs):
        self.method = method
        self.allow_empty = allow_empty
        self.parser = parser
        super(WebRequestHelper, self).__init__(**kwargs)

    def action(self, url, **kwargs):
        method = getattr(requests, self.method)
        logger.info("Sending request to %s", url)
        response = method(url, **kwargs)
        response.raise_for_status()

        if not self.allow_empty and not response.text and not self.allow_fail:
            logger.warning("Response body was empty. Retrying.")
            raise WebRequestFailed("Got empty response")
        if self.parser:
            return self.parser(response)
        return response.text


class StructureFileFinder(object):
    extension = None

    def __init__(self, config, extension=None, strict=True):
        self.config = config
        self.strict = strict
        self.location = os.path.join(self.config['locations']['fr3d_root'],
                                     'PDBFiles')
        if extension:
            self.extension = extension

        if not self.extension:
            raise ValueError("Must define an extension")

    def __call__(self, pdb):
        filename = os.path.abspath(os.path.join(self.location,
                                                pdb + '.' + self.extension))

        if not os.path.exists(filename):
            msg = "Could not find CIF file for %s.  Expected at: %s"
            raise MissingFileException(msg % (pdb, filename))
        if not self.strict:
            filename = None

        return filename


class CifFileFinder(StructureFileFinder):
    """This is a class that can find cif files. It looks into the configure
    PDBFiles directory for cif files. It will return the full path to the file,
    or raise an exception if it could not be found in strict mode. If not
    strict it returns None.
    """
    extension = 'cif'


class PickleFileFinder(StructureFileFinder):
    extension = 'pickle'
    pass


class FTPFetchHelper(RetryHelper):
    def __init__(self, uri, allow_empty=False, parser=None, **kwargs):
        self.parser = parser
        self.allow_empty = allow_empty
        logger.debug("Connecting to %s", uri)
        self.ftp = FTP(uri)
        self.ftp.login()
        logger.debug("Connected and logged in")
        super(FTPFetchHelper, self).__init__(**kwargs)

    def action(self, filename, **kwargs):
        logger.info("Attempting to get %s", filename)
        out = sio.StringIO()
        self.ftp.retrbinary("RETR %s" % filename, out.write)
        logger.info("Fetched %s", filename)
        text = out.getvalue()

        if not self.allow_empty and not text:
            logger.warning("Retrived content was empty.")
            raise EmptyResponse("Got empty content.")

        if self.parser:
            return self.parser(text)
        return text


class GzipFetchHelper(WebRequestHelper):
    def __init__(self, **kwargs):
        super(GzipFetchHelper, self).__init__(parser=self.parser)

    def parser(self, response):
        fileobj = sio.StringIO(response.content)
        unzipped = gzip.GzipFile(fileobj=fileobj)
        return unzipped.read()
