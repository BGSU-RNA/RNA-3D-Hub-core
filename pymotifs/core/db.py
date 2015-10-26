import logging
from contextlib import contextmanager

from pymotifs.core.exceptions import Skip


class Session(object):
    """A wrapper around a session maker to provide the types of logging and
    rollbacks that we desire.
    """

    def __init__(self, session_maker):
        self.logger = logging.getLogger('core.Session')
        self.maker = session_maker
        if isinstance(session_maker, Session):
            self.maker = session_maker.maker

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
        except Skip:
            session.rollback()
            raise
        except Exception:
            self.logger.error("Transaction failed. Rolling back.")
            session.rollback()
            raise
        finally:
            session.close()
