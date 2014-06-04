import os
import sys
import logging
import traceback
from contextlib import contextmanager


class MissingFileException(Exception):
    pass


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


def main(klass):
    import sys
    from models import Session
    loader = klass(Session)
    loader(sys.argv[1:])
