import os
import sys
import logging
import traceback
import itertools as it
from contextlib import contextmanager

# import requests


class MissingFileException(Exception):
    pass


def grouper(n, iterable):
    iterator = iter(iterable)
    while True:
        chunk = tuple(it.islice(iterator, n))
        if not chunk:
            return
        yield chunk


# class WebRequestHelper(object):

#     retries = 3

#     def __init__(self, method='get', retries=None):
#         self.retries = retries or self.__class__.retries
#         self.method = method

#     def __call__(self, url, **kwargs):
#         method = getattr(requests, self.method)

#         logging.info("Sending request to %s", url)
#         for index in xrange(self.retries):
#             try:
#                 response = method(url, **kwargs)
#                 response.raise_for_status()
#                 return response.text
#             except:
#                 logging.warning("Failed attempt #%s for %s", str(index), url)

#         logging.error("All attempts at fetching %s failed", url)
#         return None


class DatabaseHelper(object):

    insert_max = 1000

    def __init__(self, maker):
        self.maker = maker

    def store(self, data):
        with self.session() as session:
            for index, datum in enumerate(data):
                session.insert(datum)
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
            logging.error(traceback.format_exc(sys.exc_info()))
            logging.error("Transaction failed. Rolling back.")
            session.rollback()
            raise
        finally:
            session.close()


class CifFileFinder(object):

    def __call__(self, pdb):
        filename = os.path.join(self.config['locations']['fr3d_root'], 'FR3D',
                                'PDBFiles', pdb + '.cif')
        if not os.path.exists(filename):
            logging.warning("Could not find CIF file for %s. Expected at: %s",
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
