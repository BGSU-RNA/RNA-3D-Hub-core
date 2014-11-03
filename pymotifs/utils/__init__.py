import os
import sys
import logging
import traceback
from ftplib import FTP
import cStringIO as sio
import itertools as it
import collections as coll
from contextlib import contextmanager
import inspect
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
    iterator = iter(iterable)
    while True:
        chunk = tuple(it.islice(iterator, n))
        if not chunk:
            return
        yield chunk


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
    def __init__(self, retries=3, allow_fail=False):
        self.retries = retries
        self.allow_fail = allow_fail

    def __call__(self, *args, **kwargs):
        for index in xrange(self.retries):
            try:
                return self.action(*args, **kwargs)
            except:
                if not self.allow_fail:
                    logger.warning("Failed retry attempt #%s", str(index + 1))
                    logger.warning(traceback.format_exc(sys.exc_info()))

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

        if not self.allow_empty and not response.text:
            logger.warning("Response body was empty. Retrying.")
            raise WebRequestFailed("Got empty response")
        if self.parser:
            return self.parser(response)
        return response.text


class DatabaseHelper(object):

    insert_max = 1000

    def __init__(self, maker):
        self.maker = maker

    def store(self, data):
        with self.session() as session:
            if not isinstance(data, coll.Iterable):
                session.add(data)
            else:
                iterator = enumerate(data)
                if inspect.isgenerator(data) or \
                   isinstance(data[0], coll.Iterable) or \
                   inspect.isgenerator(data[0]):
                    iterator = enumerate(it.chain.from_iterable(data))

                for index, datum in iterator:
                    session.add(datum)
                    if index % self.insert_max == 0:
                        logger.debug("Committing a chunk of %s",
                                     self.insert_max)
                        session.commit()
                logger.debug("Final commit")

            session.commit()

    @contextmanager
    def session(self):
        session = self.maker()
        try:
            yield session
            session.commit()
        except:
            logger.error(traceback.format_exc(sys.exc_info()))
            logger.error("Transaction failed. Rolling back.")
            session.rollback()
            raise
        finally:
            session.close()


class CifFileFinder(object):

    def __init__(self, config):
        self.config = config

    def __call__(self, pdb):
        filename = os.path.join(self.config['locations']['fr3d_root'], 'FR3D',
                                'PDBFiles', pdb + '.cif')
        if not os.path.exists(filename):
            logger.warning("Could not find CIF file for %s. Expected at: %s",
                           pdb, filename)
            raise MissingFileException()
        return filename


class RNAPdbsHelper(object):
    url = 'http://www.rcsb.org/pdb/rest/search'

    def parse_all_rna_pdbs(self, raw):
        return filter(lambda x: len(x) == 4, raw.split("\n"))

    def __call__(self):
        """Get a list of all rna-containing pdb files, including hybrids. Raise
           a specific error if it fails.
        """

        logger.debug('Getting a list of all rna-containing pdbs')
        query_text = """
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>I</containsProtein>
        <containsDna>I</containsDna>
        <containsRna>Y</containsRna>
        <containsHybrid>I</containsHybrid>
        </orgPdbQuery>
        """
        post = WebRequestHelper(method='post', parser=self.parse_all_rna_pdbs)

        try:
            return post(self.url, data=query_text)
        except:
            logger.critical("Failed to get list of RNA containing PDBs")
            raise GetAllRnaPdbsError()


def row2dict(row):
    d = {}
    for column in row.__table__.columns:
        d[column.name] = str(getattr(row, column.name))
    return d


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
