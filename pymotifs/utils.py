import os
import sys
import logging
import argparse
import traceback
import itertools as it
from contextlib import contextmanager

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


def grouper(n, iterable):
    iterator = iter(iterable)
    while True:
        chunk = tuple(it.islice(iterator, n))
        if not chunk:
            return
        yield chunk


class WebRequestHelper(object):

    """A class to help with making web requests. This deals with the retrying
    and making sure the response body is not empty. If the max number of
    retries is reached we raise an exception. In addition, all steps are
    logged.
    """

    def __init__(self, allow_empty=False, method='get', retries=3,
                 parser=None):
        self.retries = retries
        self.method = method
        self.allow_empty = allow_empty
        self.parser = parser

    def __call__(self, url, **kwargs):
        method = getattr(requests, self.method)

        logger.info("Sending request to %s", url)
        for index in xrange(self.retries):
            try:
                response = method(url, **kwargs)
                response.raise_for_status()
                if not self.allow_empty and not response.text:
                    logger.warning("Response body was empty. Retrying.")
                    continue
                if self.parser:
                    return self.parser(response.text)
                return response.text
            except:
                logger.warning("Failed attempt #%s for %s", str(index + 1),
                               url)
                logger.warning(traceback.format_exc(sys.exc_info()))

        logger.error("All attempts at fetching %s failed", url)
        raise WebRequestFailed("Failed getting %s" % url)


class DatabaseHelper(object):

    insert_max = 1000

    def __init__(self, maker):
        self.maker = maker

    def store(self, data):
        with self.session() as session:
            for index, datum in enumerate(data):
                session.add(datum)
                if index % self.insert_max == 0:
                    session.commit()
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
