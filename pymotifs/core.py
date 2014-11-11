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


class Skip(Exception):
    """Base class for skipping things.
    """
    pass


class SkipPdb(Skip):
    """This is an exception that is raised to indicate we should skip
    processing this pdb. This means all values from the transformed pdb will be
    skipped as well. This can be raised during transform, has_data, and data.
    """
    pass


class SkipValue(Skip):
    """
    This is raised to indicate that we should skip the transformed value, but
    not the entire PDB. This can be raised during has_data, and data.
    """
    pass


class Matlab(object):
    """A simple wrapper around mlab. This is useful because it sets the root as
    well as calling the setup function before running matlab.
    """

    def __init__(self, root):
        self.logger = logging.getLogger('core.Matlab')
        self.logger.debug('Starting up matlab')
        os.chdir(root)
        self.mlab = mlab
        self.mlab._autosync_dirs = False
        self.mlab.setup()
        self.logger.debug('Matlab started')

    def __getattr__(self, key):
        self.logger.deug("Running %s", key)
        return getattr(self.mlab, key)


class Session(object):
    """A wrapper around a session maker to provide the types of logging and
    rollbacks that we desire.
    """

    def __init__(self, session_maker):
        self.logger = logging.getLogger('core.Session')
        self.maker = session_maker

    @contextmanager
    def __call__(self, log_exceptions=True):
        """Context handler for the session. This creates a new session and
        yields it. It will catch, log, and re raise any exceptions that occur.
        It will also commit and close all sessions.
        """
        session = self.maker()
        try:
            yield session
            session.commit()
        except Skip:
            session.rollback()
            raise
        except:
            if log_exceptions:
                self.logger.error("Transaction failed. Rolling back.")
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
            raise AttributeError("Must set update gap")

        self.session = Session(session_maker)
        self.logger = logging.getLogger(self.name)

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

    def transform(self, pdb, **kwargs):
        """This method takes the a pdb that we are given to process and
        transforms into values that can be given to `self.data`. By default it
        simply yields the given object. This is intended to be used in cases
        where we need to limit the list of pdbs to select ones. Or where we
        take list of pdbs and get a list of pairs, like with getting nt-nt
        correspondence.

        :pdb: The pdb.
        :returns: A list of things to send to data.
        """
        return [pdb]

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

        self.logger.debug("Storing data")
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
                        session.commit()

            session.commit()
            self.logger.debug("Done committing")

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

    def must_recompute(self, pdb, recalculate=False, **kwargs):
        """Detect if we have been told to recompute this stage for this pdb.
        """
        return recalculate or self.config[self.name].get('recompute')

    def should_compute(self, pdb, recalculate=False, **kwargs):
        """Determine if we should recompute the data for this loader and pdb.
        This is true if we are told to recompute, if we do not have data for
        this pdb or it has been long enough since the last update.
        """
        if self.must_recompute(pdb, recalculate=recalculate, **kwargs):
            self.logger.debug("Performing a forced recompute")
            return True

        has_data = self.has_data(pdb)
        if not has_data:
            self.logger.debug("Missing data for %s will compute", pdb)
            return True

        too_long = self.been_long_enough(pdb)
        if too_long:
            self.logger.debug("Time gap for %s too large, recomputing", pdb)
            return True

        self.logger.debug("No reason to recompute %s", pdb)
        return False

    def mark_analyzed(self, pdb):
        """Mark that we have finished computing the results for the given pdb.
        """
        with self.session() as session:
            status = PdbAnalysisStatus(id=pdb, step=self.step(),
                                       time=datetime.datetime.now())
            session.merge(status)
            session.commit()
            self.logger.info('Updated %s status for pdb %s', self.name, pdb)

    def __call__(self, pdbs, **kwargs):
        """Load all data for the list of pdbs. This will load things as needed
        into the database. It checks if we should recompute or if we are forced
        to and then stores the computed data.

        :pdbs: A list of pdbs to process.
        :kwargs: Keyword arguments. These are passed along to other methods.
        """

        if not pdbs:
            self.logger.critical("Must give pdbs")
            raise InvalidState("No pdbs given")

        for index, pdb in enumerate(pdbs):
            self.logger.info("Starting %s: %s/%s", pdb, index + 1, len(pdbs))
            pdb = pdb.upper()

            try:
                iterable = self.transform(pdb, **kwargs)
            except SkipPdb as err:
                self.logger.warn("Skipping pdb %s. Reason: %s", pdb, str(err))
                continue
            except:
                self.logger.error("Transforming %s failed, skipping", pdb)
                self.logger.error(traceback.format_exc(sys.exc_info()))
                continue

            if not iterable:
                self.logger.info("Nothing from transfrom for %s", pdb)
                continue

            for trans_index, transformed in enumerate(iterable):
                self.logger.info("Processing %s (%s) %s: %s/%s",
                                 transformed, pdb, index, trans_index + 1,
                                 len(iterable))

                try:
                    if not self.should_compute(transformed, **kwargs):
                        self.logger.debug("No need to compute %s", transformed)
                        continue
                except SkipPdb as err:
                    self.logger.warn("Skipping pdb %s Reason %s",
                                     pdb, str(err))
                    break
                except SkipValue as err:
                    self.logger.warn("Skipping %s Reason %s",
                                     transformed, str(err))
                    continue

                if self.must_recompute(transformed, **kwargs):
                    self.logger.debug("Removing old data for %s", transformed)
                    self.remove(transformed)

                try:
                    data = self.data(transformed, **kwargs)
                except SkipPdb as err:
                    self.logger.warn("Skipping pdb: %s Reason %s",
                                     pdb, str(err))
                    break
                except SkipValue as err:
                    self.logger.warn("Skipping %s Reason %s",
                                     transformed, str(err))
                    continue
                except:
                    if self.stop_on_failure:
                        self.logger.error("Error raised when getting data")
                        self.logger.error(traceback.format_exc(sys.exc_info()))
                        raise StageFailed(self.name)

                    self.logger.warning("Error raised when getting data")
                    self.logger.warning(traceback.format_exc(sys.exc_info()))
                    continue

                if not self.allow_no_data and not data:
                    self.logger.error("No data produced")
                    raise InvalidState("Missing data")
                elif not data:
                    self.logger.warning("No data produced")
                    continue

                try:
                    self.store(data)
                except:
                    self.remove(transformed)
                    if self.stop_on_failure:
                        self.logger.error("Error raised when storing data")
                        self.logger.error(traceback.format_exc(sys.exc_info()))
                        raise StageFailed(self.name)

                    self.logger.warning("Error raised when storing data")
                    self.logger.warning(traceback.format_exc(sys.exc_info()))
                    continue


class MultiLoader(object):
    stages = []

    def __init__(self, *args):
        self.steps = []
        for stage in self.stages:
            self.steps.append(stage(*args))

    def __call__(self, pdbs, **kwargs):
        for step in self.steps:
            step(pdbs, **kwargs)
