"""This modules contains the classes for building stages of the pipeline.
"""

import os
import abc
import pickle
import datetime

from fr3d.data import Structure
from fr3d.cif.reader import Cif
from fr3d.cif.reader import ComplexOperatorException

from pymotifs.core import base
from pymotifs.core.exceptions import Skip
from pymotifs.core.exceptions import StageFailed
from pymotifs.core.exceptions import InvalidState
from pymotifs import utils as ut
from pymotifs import models as mod
from pymotifs.core import savers

# This is a very large virus file that should be skipped. Add other files as
# necessary
SKIP = set(['4V3P', '4V4G'])


class Stage(base.Base):
    """This is a base class for both loaders and exporters to inherit from. It
    contains the functionality common to all things that are part of our
    pipeline.
    """

    """ Maximum length of time between updates. False for forever. """
    update_gap = None

    """ What stages this stage depends upon. """
    dependencies = set()

    """Flag if we should mark stuff as processed."""
    mark = True

    """If we should skip complex operators"""
    skip_complex = True

    """Collection of ids to skip"""
    skip = []

    """Class to save with"""
    saver = None

    def __init__(self, *args, **kwargs):
        """Build a new Stage.

        :config: The config object to build with.
        :session_maker: The Session object to handle database connections.
        """
        super(Stage, self).__init__(*args, **kwargs)
        self._cif = ut.CifFileFinder(self.config)
        self.skip = set(SKIP)
        self.skip.update(self.__class__.skip)
        self.skip.update(kwargs.get('skip_pdbs', []))

    @abc.abstractmethod
    def is_missing(self, entry, **kwargs):
        """Determine if we do not have any data. If we have no data then we
        will recompute.

        :entry: The thing to check for.
        :kwargs: Generic keyword arguments
        :returns: True or False
        """
        pass

    @abc.abstractmethod
    def process(self, entry, **kwargs):
        """Process this entry. In the case of loaders this will parse the data
        and put it into the database, exporters may go to the database and then
        generate the file.

        :entry: The entry to process.
        :kwargs: Generic keyword arguments.
        :returns: Nothing and is ignored.
        """
        pass

    def remove(self, entry, **kwargs):
        pass

    def cif(self, pdb):
        """A method to load the cif file for a given pdb id.

        :pdb: PDB id to parse.
        :returns: A parsed cif file.
        """
        if isinstance(pdb, Cif):
            return pdb

        try:
            with open(self._cif(pdb), 'rb') as raw:
                return Cif(raw)
        except ComplexOperatorException as err:
            if self.skip_complex:
                self.logger.warning("Got a complex operator for %s, skipping",
                                    pdb)
                raise Skip("Complex operator must be skipped")
            raise err

    def structure(self, pdb):
        """A method to load the cif file and get the structure for the given

        :pdb: The pdb id to get the structure for.
        :returns: The FR3D structure for the given PDB.
        """
        if isinstance(pdb, Structure):
            return pdb

        return self.cif(pdb).structure()

    def cache_filename(self, name):
        cache_dir = self.config['locations']['cache']
        if not os.path.isdir(cache_dir):
            os.mkdir(cache_dir)

        return os.path.join(cache_dir, name + '.pickle')

    def cache(self, name, data):
        """Cache some data under a name.

        :name: The name to cache under.
        :data: The data to cache.
        """

        filename = self.cache_filename(name)
        with open(filename, 'wb') as raw:
            pickle.dump(data, raw)

    def cached(self, name, remove=False):
        """Load some cached data.

        :name: Name of the cached data.
        :returns: The data or None if no cached data of that name exists
        """

        filename = self.cache_filename(name)
        data = None
        if os.path.exists(filename):
            with open(filename, 'rb') as raw:
                data = pickle.load(raw)

        if remove:
            if not os.path.exists(filename):
                self.logger.warning("Attempt to remove nonexisting cache %s",
                                    name)
            else:
                os.remove(filename)

        return data

    def must_recompute(self, pdb, recalculate=False, **kwargs):
        """Detect if we have been told to recompute this stage for this pdb.
        """
        if recalculate is True:
            return True
        if isinstance(recalculate, (set, list, tuple)):
            try:
                if self.name in recalculate:
                    return True
            except:
                pass
        return bool(self.config[self.name].get('recompute'))

    def been_long_enough(self, pdb, ignore_time=False, **kwargs):
        """Determine if it has been long enough to recompute the data for the
        given pdb. This uses the udpate_gap property which tells how long to
        wait between updates. If that is False then we never update based upon
        time.
        """
        if not self.update_gap or ignore_time:
            return False

        with self.session() as session:
            current = session.query(mod.PdbAnalysisStatus).\
                filter_by(pdb_id=pdb, stage=self.name).\
                first()
            if not current:
                return True
            current = current.time
        # If this has been marked as done in the far future do it anyway. That
        # is a silly thing to do
        diff = abs(datetime.datetime.now() - current)
        return diff > self.update_gap

    def should_process(self, entry, **kwargs):
        """Determine if we should process this entry. This is true if we are
        told to recompute, if we do not have data for this pdb or it has been
        long enough since the last update.

        :entry: The entry to check.
        :kwargs: Some keyword arguments for determining if we should process
        :returns: True or False
        """
        try:
            if entry in self.skip:
                raise Skip("Forced skip of %s", entry)
        except Skip as err:
            raise err
        except Exception:
            pass

        if self.must_recompute(entry, **kwargs):
            self.logger.debug("Performing a forced recompute")
            return True

        if self.been_long_enough(entry, **kwargs):
            self.logger.debug("Time gap for %s too large, recomputing", entry)
            return True

        if self.is_missing(entry, **kwargs):
            self.logger.debug("Missing data from %s. Will compute", entry)
            return True
        return False

    def to_process(self, pdbs, **kwargs):
        """Compute the things to process. For things that work with PDBs the
        default one will work well. For things that work with groups this could
        simply ignore the given pdbs and return the names of the groups to work
        with.

        :pdbs: Input pdbs
        :kwargs: Generic keyword arguments.
        :returns: The stuff to process.
        """
        return [str(p).upper() for p in pdbs if p.upper() not in self.skip]

    def mark_processed(self, pdb, dry_run=False, **kwargs):
        """Mark that we have finished computing the results for the given pdb.

        :pdb: The pdb to mark done.
        """

        if dry_run:
            self.logger.debug("Marking %s as done", pdb)
        else:
            with self.session() as session:
                status = mod.PdbAnalysisStatus(pdb_id=pdb, stage=self.name,
                                               time=datetime.datetime.now())
                session.merge(status)
        self.logger.info('Updated %s status for pdb %s', self.name, pdb)

    def __call__(self, given, **kwargs):
        """Process all given inputs. This will first transform all inputs with
        the `to_process` method. If there are no entries then a critical
        exception is raised. We then use `should_process` to determine if we
        should process each entry. If this returns
        true then we call `process`. Once done we call `mark_processed`.

        :given: A list of pdbs to process.
        :kwargs: Keyword arguments passed on to various methods.
        :returns: Nothing
        """

        entries = None
        try:
            entries = self.to_process(given, **kwargs)
        except Skip as err:
            self.logger.warn("Skipping this stage. Reason %s", str(err))
            return []

        if not entries:
            self.logger.critical("Nothing to process")
            raise InvalidState("Nothing to process")

        failed = []
        processed = []
        for index, entry in enumerate(entries):
            self.logger.info("Processing %s: %s/%s", entry, index + 1,
                             len(entries))

            try:
                if not self.should_process(entry, **kwargs):
                    self.logger.debug("No need to process %s", entry)
                    continue
                self.process(entry, **kwargs)

            except Skip as err:
                self.logger.warn("Skipping entry %s. Reason %s",
                                 str(entry), str(err))
                continue

            except Exception as err:
                self.logger.error("Error raised in processing of %s", entry)
                self.logger.exception(err)

                try:
                    self.remove(entry, **kwargs)
                except Exception as err:
                    raise InvalidState("Could not cleanup failed data %s",
                                       entry)
                else:
                    failed.append(entry)
                    continue

            if self.mark:
                self.mark_processed(entry, **kwargs)
            processed.append(entry)

        if failed:
            ids = ', '.join(str(f) for f in failed)
            raise StageFailed("Pipeline failed on these inputs %s", ids)

        return processed


class Loader(Stage):
    """An abstract baseclass for all things that load data into our database.
    This provides a constituent interface for all loaders to use. This extends
    Stage by making the process method
    """

    __metaclass__ = abc.ABCMeta

    """ Max number of things to insert at once. """
    insert_max = 1000

    """ A flag to indicate it is ok to produce no data. """
    allow_no_data = False

    """ A flag to indicate if we should use sessions .merge instead of .add """
    merge_data = False

    """Sqlalchmey model to save to"""
    table = None

    """Class to use for saving"""
    saver = savers.DatabaseSaver

    @abc.abstractmethod
    def data(self, pdb, **kwargs):
        """Compute the data for the given cif file.
        """
        pass

    @abc.abstractmethod
    def has_data(self, pdb, **kwargs):
        """Check if we have already stored data for this pdb file in the
        database. This is used to determine if we should attempt to compute new
        data.
        """
        pass

    @abc.abstractmethod
    def remove(self, pdb, **kwargs):
        """Remove any old data for this pdb. This is used both when we have
        failed to store all data as well as when we are overwriting existing
        data. This should never raise anything unless something is drastically
        wrong.
        """
        pass

    def is_missing(self, entry, **kwargs):
        """Determine if the data is missing by using the has_data method.

        :entry: The entry to check for.
        :kwargs: Keyword arguments
        :returns: A boolean if the requested data is missing or not.
        """
        return not self.has_data(entry, **kwargs)

    def store(self, pdb, data, **kwargs):
        """Store the given data. The data is written in chunks of
        self.insert_max at a time. The data can be a list or a nested set of
        iterables. If dry_run is true then this will not actually store
        anything, but instead will log the attempt. If this loader has
        'merge_data' set to True then this will merge instead of adding data.

        :data: The data to store. May be a list, an iterable nested one level
        deep or a single object to store.
        :dry_run: A flag to indicate if this should perform a dry run.
        :kwargs: Keyword arguments.
        """
        saver = self.saver(self.config, self.session, stage=self)
        saver(pdb, data, **kwargs)

    def process(self, entry, **kwargs):
        """Get the data for a particular entry. This will get the data and then
        store it. It will remove data as needed and makes sure that data is
        produced if required.

        :entry: The entry to process.
        :kwargs: Generic keyword arguments to be passed along to other methods.
        """

        if self.must_recompute(entry, **kwargs):
            if kwargs.get('dry_run'):
                self.logger.debug("Skipping removal in dry run")
            else:
                self.logger.debug("Removing old data for %s", entry)
                self.remove(entry)

        data = self.data(entry, **kwargs)

        if not data:
            if not self.allow_no_data:
                self.logger.error("No data produced")
                raise InvalidState("Stage %s produced no data processing %s" %
                                   (self.name, entry))
            else:
                self.logger.warning("No data produced")
                return

        self.store(entry, data, **kwargs)


class SimpleLoader(Loader):
    """A SimpleLoader is a subclass of Loader that has a default implementation
    of has_data and remove. These depend on the abstract method query which
    generates the query to use for these things. Basically this is what to use
    if we are simply adding things to a table in the database.
    """

    __metaclass__ = abc.ABCMeta

    def has_data(self, args, **kwargs):
        """Check if we already have data.
        """
        with self.session() as session:
            return bool(self.query(session, args).limit(1).count())

    def remove(self, args, **kwargs):
        """This will delete all entries for the given arguments
        """

        self.logger.info("Removing data for %s", str(args))
        if kwargs.get('dry_run'):
            return True

        with self.session() as session:
            query = self.query(session, args)
            if not query.count():
                self.logger.info("Nothing to delete for %s", str(args))
                return True

            for row in query:
                session.delete(row)

    @abc.abstractmethod
    def query(self, session, args):
        """
        A method to generate the query that can be used to access data for this
        loader. The resutling query is used in remove and has_data.

        :session: The session object to use.
        :*args: Arguments from process.
        """
        pass


class MassLoader(Loader):
    """A MassLoader is a Loader that works on collections of PDB files. For
    example, when getting all PDB info we do that for all PDB files at once in
    a single request. Here we do many of the normal things that a stage does
    but we simply do it on all pdbs at once.
    """

    __metaclass__ = abc.ABCMeta

    def been_long_enough(self, pdbs, **kwargs):
        """Determine if it has been long enough to recompute the data for the
        given pdbs. This requires that it has been long enough for at least 1
        structure.
        """
        parent = super(MassLoader, self).been_long_enough
        return any(parent(pdb, **kwargs) for pdb in pdbs)

    def to_process(self, pdbs, **kwargs):
        return [tuple(super(MassLoader, self).to_process(pdbs))]

    def should_process(self, pdbs, **kwargs):
        """We check if we have data for all pdbs. If we are missing 1 then we
        will recompute.
        """
        parent = super(MassLoader, self).should_process
        return any(parent(pdb, **kwargs) for pdb in pdbs)

    def mark_processed(self, pdbs, **kwargs):
        for pdb in pdbs:
            super(MassLoader, self).mark_processed(pdb, **kwargs)

    def remove(self, pdbs, **kwargs):
        self.logger.debug("Remove does nothing in MassLoaders")

    @abc.abstractmethod
    def data(self, pdbs, **kwargs):
        pass


class StageContainer(Stage):
    """This acts as a simple way to aggregate a bunch of loaders into one
    stage. It is really useful sometimes to simply run say all unit loaders
    without having to do each one individually or know which ones depend on
    each other. The loader itself does nothing but provide a way to run all
    other loaders this depends on. All stages which inherit from this will do
    nothing by themselves other than to run other stages. Do not try to add
    behavior to these loaders.
    """

    """The list of stages that are children of this loader"""
    stages = []


class Exporter(Loader):
    """A class that saves to CSV files.
    """
    __metaclass__ = abc.ABCMeta

    """List of headers to save."""
    headers = []

    """Class to use for writing to csv files."""
    saver = savers.CsvSaver

    @abc.abstractmethod
    def filename(self, entry, **kwargs):
        """Compute the filename for the given entry.

        :entry: The entry to write out.
        """
        pass

    def to_process(self, pdbs, **kwargs):
        """If we are given the all flag then we should get all known pdb ids,
        otherwise just use the ones we are given. This will turn all the given
        pdbs into a list of lists, so they are all processed in one go.

        :pdbs: The list of pdbs to get interactions for.
        :kwargs: All keyword arguments.
        :returns: The list of pdbs to get interactions for.
        """

        if not kwargs.get('all'):
            return [tuple(pdbs)]

        with self.session() as session:
            query = session.query(mod.PdbInfo.pdb_id).distinct()
            return [tuple([result.pdb_id for result in query])]

    def has_data(self, *args, **kwargs):
        """We always recompute when exporting.
        """
        return False

    def remove(self, *args, **kwargs):
        """We never remove exported files automatically.
        """
        self.logger.info("No automatic removal in exporters")
