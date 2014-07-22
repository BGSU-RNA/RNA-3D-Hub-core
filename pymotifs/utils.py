import os
import sys
import logging
import traceback
import itertools as it
from contextlib import contextmanager

import requests

# logger = logging.getLogger('utils')
logger = logging


class MissingFileException(Exception):

    """This is raised when we can't find a file. For example a cif file for a
    pdb does not exist.
    """
    pass


class WebRequestFailed(Exception):

    """Raised when we have failed all attempts at getting a url.
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

    def __init__(self, allow_empty=False, method='get', retries=3):
        self.retries = retries
        self.method = method
        self.allow_empty = allow_empty

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
                return response.text
            except:
                logger.warning("Failed attempt #%s for %s", str(index), url)

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
    import sys
    from models import Session
    loader = klass(Session)
    loader(sys.argv[1:])
