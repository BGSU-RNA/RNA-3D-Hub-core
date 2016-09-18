"""This module contains classes for saving data. We have to save data into two
general places, files and the database. This was created to provide a constient
interface for saving to the stages. These classes also deal with error handling
for when things can't get written, or writing fails for some reason.

In addition, these contain logic for doing a 'dry run', that is running the
pipeline without actually saving anything which can be useful for debugging
specific stages.
"""

import os
import abc
import csv
import gzip
import shutil
import collections as coll
from contextlib import contextmanager

from pymotifs import utils as ut

from pymotifs.core.base import Base
from pymotifs.core.exceptions import InvalidState


class SaveFailed(Exception):
    """This is raised when there is an issue with saving. For example if we
    write to a file but nothing is created, this raises an exception.
    """
    pass


class Saver(Base):
    """This is a base class for saving. This cannot be used directly, it is
    just a base class for all other savers.

    Attributes
    ----------
    stage : pymotifs.core.stages.Stage
        The stage this saver is a part of
    allow_no_data : bool
        A flag to control if we should fail if no data is written
    merge : bool
        If we should merge instead of saveing normally.
    insert_max : int
        The maximum number of entries to write at once.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, config, session, stage=None, **kwargs):
        """Create a new saver.

        Parameters
        ----------
        config : dict
            The configuration dict
        session : pymotifs.core.db.Session
            The session wrapper
        stage : pymotifs.core.stages.Stage
            The stage this saver is a part of
        """

        super(Saver, self).__init__(config, session)
        self.stage = stage or object()
        self.allow_no_data = getattr(stage, 'allow_no_data', False)
        self.merge = getattr(stage, 'merge_data', False)
        self.insert_max = getattr(stage, 'insert_max', 1000)

    @abc.abstractmethod
    def writer(self, pdb, **kwargs):
        """Create the writer to use. Classes must implement this and it must
        be a context handler. It should produce a function which may be called
        to save data.
        """
        pass

    @contextmanager
    def _writer(self, pdb, **kwargs):
        """Create the writer to use. If given the keyword argument 'dry_run' and
        it is set to try then we will yield a function which simply logs the
        data given. Otherwise this calls the .writer method, which must be
        implemented by inheriting classes to create the writer.
        """

        if kwargs.get('dry_run'):
            yield lambda data: self.logger.debug("Saving %s", data)
        else:
            with self.writer(pdb, **kwargs) as writer:
                yield writer

    def __call__(self, pdb, data, **kwargs):
        """Save data for the given pdb.

        Parameters
        ----------
        pdb : str
            The pdb this a part of
        data : list or dict
            The data to save
        """

        to_save = data
        if not isinstance(to_save, coll.Iterable) or isinstance(to_save, dict):
            to_save = [data]

        saved = False
        for index, chunk in enumerate(ut.grouper(self.insert_max, to_save)):
            chunk = list(chunk)
            kwargs['index'] = index
            with self._writer(pdb, **kwargs) as writer:
                for entry in chunk:
                    writer(entry)
                    saved = True

        if not saved:
            if not self.allow_no_data:
                raise InvalidState("No data saved")
            self.logger.warning("No data saved")


class DatabaseSaver(Saver):
    """A saver that writes to the database. This uses the Session wrapper to
    write to the database. This can deal with either inserting or updating.
    This is controlled by setting merge to True on the Stage this belongs to
    for updating.

    If `table` is set on the stage then it is possible to save a dictonary by
    having it automatically converted to the correct object. This can be very
    useful for logging and debugging as the data remains in a nice printable
    form until the last possible moment.
    """

    def __init__(self, *args, **kwargs):
        super(DatabaseSaver, self).__init__(*args, **kwargs)
        self.table = getattr(self.stage, 'table', None)

    def to_savable(self, data):
        """Turn data into a table object if needed.
        """
        if isinstance(data, dict) and self.table:
            return self.table(**data)
        return data

    @contextmanager
    def writer(self, *args, **kwargs):
        with self.session() as session:
            fn = session.add
            if self.merge:
                fn = session.merge
            yield lambda data: fn(self.to_savable(data))


class FileHandleSaver(Saver):
    """A saver that produces a file handle as a writer. This is intended to be
    inherited from for creating new savers. This can't be used directly. It
    also has the ability to compress files after writing if needed. This will
    also check that the file is created and is not empty after writing. If that
    happens then this will raise an exception if allow_no_data is False.
    """

    def __init__(self, *args, **kwargs):
        super(FileHandleSaver, self).__init__(*args, **kwargs)
        self.compressed = getattr(self.stage, 'compressed', False)

    def filename(self, pdb, **kwargs):
        """This will use the stage's filename method for computing the filename
        to save to.
        """
        return self.stage.filename(pdb, **kwargs)

    def mode(self, pdb, **kwargs):
        if not self.merge and not kwargs.get('index'):
            return 'wb'
        return 'ab'

    @contextmanager
    def writer(self, pdb, **kwargs):
        """Creates a new writer to save to.
        """
        filename = self.filename(pdb, **kwargs)
        mode = self.mode(pdb, **kwargs)
        with open(filename, mode) as raw:
            yield raw

        if not os.path.isfile(filename):
            if not self.allow_no_data:
                raise SaveFailed("No file created")
            self.logger.warn("%s not created", filename)

        if not os.path.getsize(filename):
            if not self.allow_no_data:
                raise SaveFailed("Nothing written")
            self.logger.warn("Nothing written to %s", filename)

    def compress(self, pdb, data, **kwargs):
        """Compress the file for the given pdb and data.
        """

        filename = self.filename(pdb, **kwargs)
        self.logger.debug('Compressing %s', filename)
        if kwargs.get('dry_run'):
            return True

        temp_output_file = os.path.join(self.config['fr3d_root'],
                                        self.__name__)
        handle = gzip.open(temp_output_file, 'wb')
        with open(filename, 'rb') as raw:
            handle.writelines(raw)
            handle.close()
        shutil.move(temp_output_file, filename)

    def __call__(self, *args, **kwargs):
        super(FileHandleSaver, self).__call__(*args, **kwargs)
        if self.compressed:
            self.compress(*args, **kwargs)


class CsvSaver(FileHandleSaver):
    """A Saver for writing CSV files. This will save a csv file by writing the
    given dictonary using the configured headers. It will also write a header
    line. Headers are configured with the `headers` property on the stage this
    belongs to.
    """

    def __init__(self, *args, **kwargs):
        super(CsvSaver, self).__init__(*args, **kwargs)
        self.headers = getattr(self.stage, 'headers')
        self._headers_written = False

    @contextmanager
    def writer(self, pdb, **kwargs):
        """This will create a csv writer with the standard options we use.
        """
        with super(CsvSaver, self).writer(pdb, **kwargs) as handle:
            writer = csv.DictWriter(handle, self.headers, quotechar='"',
                                    quoting=csv.QUOTE_ALL)

            if self._headers_written:
                writer.writerow(dict(zip(self.headers, self.headers)))

            yield writer.writerow
