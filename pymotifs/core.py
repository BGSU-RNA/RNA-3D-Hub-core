import os
import abc
import sys
import logging
import inspect
import datetime
import traceback
import itertools as it
import collections as coll
from contextlib import contextmanager

try:
    from mlabwrap import mlab
except:
    pass

from models import PdbAnalysisStatus

logger = logging.getLogger(__name__)


class StageFailed(Exception):
    """This is raised when one stage of the pipeline fails.
    """
    pass


class InvalidState(Exception):
    """This is an exception meant to be used when we have entered into some
    sort of invalid state in a stage. For example, we require that chain breaks
    be loaded before we do nt-nt correspondences but if the programs are run
    out of order this will be raised.
    """
    pass


class Matlab(object):
    """A simple wrapper around mlab. This is useful because it sets the root as
    well as calling the setup function before running matlab.
    """

    def __init__(self, root):
        logger.debug('Starting up matlab')
        os.chdir(root)
        self.mlab = mlab
        self.mlab._autosync_dirs = False
        self.mlab.setup()
        logger.debug('Matlab started')

    def __getattr__(self, key):
        return getattr(self.mlab, key)


class Session(object):
    """A wrapper around a session maker to provide the types of logging and
    rollbacks that we desire.
    """

    def __init__(self, session_maker):
        self.maker = session_maker

    @contextmanager
    def __call__(self):
        """Context handler for the session. This creates a new session and
        yields it. It will catch, log, and re raise any exceptions that occur.
        It will also commit and close all sessions.
        """
        session = self.maker()
        try:
            yield session
            session.commit()
        except:
            logger.error("Transaction failed. Rolling back.")
            logger.error(traceback.format_exc(sys.exc_info()))
            session.rollback()
            raise
        finally:
            session.close()


class Loader(object):
    """An abstract baseclass for all loaders. This provides a constient
    interface for all loaders to use.
    """

    __metaclass__ = abc.ABCMeta

    """ The name of the loader. """
    name = None

    """ Maximum length of time between updates. False for forever. """
    update_gap = None

    """ If we should stop the whole loader if one part fails. """
    stop_on_failure = True

    """ Max number of things to insert at once. """
    insert_max = 1000

    """ A flag to indicate it is ok to produce no data. """
    allow_no_data = False

    def __init__(self, config, session_maker):
        self.config = coll.defaultdict(dict)
        self.config.update(config)
        if not self.name:
            raise AttributeError("Must set name")

        if self.update_gap is None:
            raise AttributeError("Must set update gap for %s" % self.name)

        self.session = Session(session_maker)

    @abc.abstractmethod
    def data(self, pdb, **kwargs):
        """Compute the data for the given pdb file.
        """
        pass

    @abc.abstractmethod
    def has_data(self, pdb):
        """Check if we have already stored data for this pdb file in the
        database. This is used to determine if we should attempt to compute new
        data.
        """
        pass

    @abc.abstractmethod
    def remove(self, pdb):
        """Remove any old data for this pdb. This is used both when we have
        failed to store all data as well as when we are overwriting existing
        data. This should never raise anything unless something is drastically
        wrong.
        """
        pass

    def transform(self, pdbs):
        """This method takes the list of pdbs that we are given to process and
        transfrorms into values that can be given to data. By default it simply
        returns the given list. This is intended to be used in cases where we
        need to limit the list of pdbs to select ones. Or where we take list of
        pdbs and get a list of pairs, like with getting nt-nt correspondence.

        :pdbs: The list of pdbs.
        :returns: A list of things to send to data.
        """
        return pdbs

    def step(self):
        """Gets the name of this step. Basically it is used to normalize
        self.name.
        """
        return self.name.lower()

    def store(self, data):
        """Store the given data. The data is written in chunks of
        self.insert_max at a time. The data can be a list or a nested set of
        iterables.
        """

        logger.debug("Storing data for %s", self.name)
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

    def been_long_enough(self, pdb):
        """Determine if it has been long enough to recompute the data for the
        given pdb. This uses the udpate_gap property which tells how long to
        wait between updates. If that is False then we never update based upon
        time.
        """
        if not self.update_gap:
            return False

        with self.session() as session:
            current = session.query(PdbAnalysisStatus).\
                filter_by(id=pdb, step=self.step()).\
                first()
            if not current:
                return True
            current = current.time
        # If this has been marked as done in the far future do it anyway. That
        # is a silly thing to do
        diff = abs(datetime.datetime.now() - current)
        return diff > self.update_gap()

    def should_compute(self, pdb, recalculate=False, **kwargs):
        """Determine if we should recompute the data for this loader and pdb.
        This is true if we are told to recompute, if we do not have data for
        this pdb or it has been long enough since the last update.
        """
        must_recompute = recalculate or self.config[self.name].get('recompute')
        if must_recompute:
            logger.debug("Given recompute for %s", pdb)
            return True

        has_data = self.has_data(pdb)
        if not has_data:
            logger.debug("Missing data for %s will recompute", pdb)
            return True

        too_long = self.been_long_enough(pdb)
        if too_long:
            logger.debug("Time gap for %s too large, recomputing", pdb)
            return True

        logger.debug("No reason to recompute %s", pdb)
        return False

    def mark_analyzed(self, pdb):
        """Mark that we have finished computing the results for the given pdb.
        """
        with self.session() as session:
            status = PdbAnalysisStatus(id=pdb, step=self.step(),
                                       time=datetime.datetime.now())
            session.merge(status)
            session.commit()
            logging.info('Updated %s status for pdb %s', self.name, pdb)

    def __call__(self, pdbs, **kwargs):
        """Load all data for the list of pdbs. This will load things as needed
        into the database. It checks if we should recompute or if we are forced
        to and then stores the computed data.
        """

        if not pdbs:
            logging.critical("Must give pdbs to %s", self.name)
            raise Exception("Not pdbs given")

        transformed = self.transform(pdbs)
        if not transformed:
            logging.info("Nothing returned from transformation")
            return None

        failed_count = 0
        for pdb in transformed:
            pdb = pdb.upper()
            logger.info("Loader %s is processing %s", self.name, pdb)

            if not self.should_compute(pdb, **kwargs):
                logger.debug("Skipping computing %s for %s", self.name, pdb)
                continue

            try:
                data = self.data(pdb, **kwargs)
            except:
                logger.error("Error raised when getting data for %s",
                             self.name)
                logger.error(traceback.format_exc(sys.exc_info()))
                self.remove(pdb)
                if self.stop_on_failure:
                    raise StageFailed(self.name)
                failed_count += 1
                continue

            if not self.allow_no_data and not data:
                logger.error("No data produced for %s", self.name)
                raise InvalidState("Missing data")
            elif not data:
                logger.warning("No data produced for %s", self.name)
                failed_count += 1
                continue

            try:
                self.remove(pdb)
                self.store(data)
            except:
                logger.error("Error raised when storing data for %s",
                             self.name)
                logger.error(traceback.format_exc(sys.exc_info()))
                self.remove(pdb)
                if self.stop_on_failure:
                    raise StageFailed(self.name)
                failed_count += 1
                continue

            # try:
            #     self.mark_analyzed(pdb)
            # except:
            #     logger.error("Could not mark %s as done on %s", self.name, pdb)
            #     logger.error(traceback.format_exc(sys.exc_info()))
            #     self.remove(pdb)
            #     if self.stop_on_failure:
            #         raise StageFailed(self.name)
            #     failed_count += 1

        logger.info("%s out of %s pdbs failed", failed_count, len(pdbs))
        if failed_count == len(pdbs):
            logger.error("All pdbs failed")


class MultiLoader(object):
    stages = []

    def __init__(self, *args):
        self.steps = []
        for stage in self.stages:
            self.steps.append(stage(*args))

    def __call__(self, pdbs, **kwargs):
        for step in self.steps:
            step(pdbs, **kwargs)
