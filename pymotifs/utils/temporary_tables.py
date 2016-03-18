"""This contains a class for working with temporary tables. It wraps up the
creation, storing and clean up of a temporary table.
"""

from contextlib import contextmanager

from pymotifs import core
from pymotifs import models as mod


class Temporary(core.Base):
    """A class to wrap working with a temporary table.
    """

    def __init__(self, table, *args, **kwargs):
        self.table = table
        self.saver = core.savers.DatabaseSaver(*args, **kwargs)
        self.saver.table = table
        super(Temporary, self).__init__(*args, **kwargs)

    @contextmanager
    def __call__(self, data=None):
        """A context manager for working with temp tables. If given data it
        will store the data in the temp table. This will produce a session
        which has the temp table and once the session ends it will delete the
        table.

        :param data: The data to store.
        """

        engine = mod.metadata.bind
        table = self.table.__table__
        try:
            table.create(engine)
            if data:
                self.saver('temp', data)
            with self.session() as session:
                yield session
        except Exception as err:
            self.logger.error("Could not work with temp table")
            self.logger.exception(err)
            raise err
        finally:
            try:
                table.drop(engine)
            except Exception as err:
                self.logger.error("Failed dropping temp")
                self.logger.exception(err)
