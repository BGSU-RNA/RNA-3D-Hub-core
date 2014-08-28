import os
import sys
import logging
import argparse
import traceback
from ftplib import FTP
import cStringIO as sio
import itertools as it
import collections as coll
from contextlib import contextmanager
import inspect

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


def append_libs():
    sys.path.append(os.path.join(os.path.dirname(__file__), "rnastructure"))


def grouper(n, iterable):
    iterator = iter(iterable)
    while True:
        chunk = tuple(it.islice(iterator, n))
        if not chunk:
            return
        yield chunk


class RetryHelper(object):
    def __init__(self, retries=3):
        self.retries = retries

    def __call__(self, *args, **kwargs):
        for index in xrange(self.retries):
            try:
                return self.action(*args, **kwargs)
            except:
                logger.warning("Failed retry attempt #%s", str(index + 1))
                logger.warning(traceback.format_exc(sys.exc_info()))

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
            return self.parser(response.text)
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


def row2dict(row):
    d = {}
    for column in row.__table__.columns:
        d[column.name] = str(getattr(row, column.name))
    return d


def main(klass):
    from models import Session

    parser = argparse.ArgumentParser(description="Run %s" % klass.__name__)
    parser.add_argument('pdbs', metavar='P', nargs='+',
                        help="PDBs to use")
    parser.add_argument('--all', dest='_all', default=False,
                        help="Use all RNA containing PDBS")
    parser.add_argument('--recalculate', action='store_true',
                        help="Force all data to be recalculated")
    parser.add_argument('--log-file', dest='_log_file', default='',
                        help="Log file to use")
    parser.add_argument('--log-level', dest='_log_level', default='debug',
                        choices=['debug', 'info', 'warning', 'error'],
                        help="Logging level to use")

    args = parser.parse_args()

    pdbs = args.pdbs
    if args._all:
        from PdbInfoLoader import PdbInfoLoader
        P = PdbInfoLoader()
        P.get_all_rna_pdbs()
        pdbs = P.pdbs

    log_args = {
        'level': getattr(logging, args._log_level.upper())
    }
    if args._log_file:
        log_args['filename'] = args._log_file

    logging.basicConfig(**log_args)

    kwargs = {}
    for arg, value in vars(args).items():
        if arg != 'pdbs' and arg[0] != '_':
            kwargs[arg] = value

    obj = klass(Session)
    obj(pdbs, **kwargs)


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
