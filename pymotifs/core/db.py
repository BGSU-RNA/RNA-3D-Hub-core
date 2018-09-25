"""This contains core classes that help with interaction with the database.
"""

import logging
from contextlib import contextmanager

from pymotifs.core.exceptions import Skip


class Session(object):
    """A wrapper around a session maker to provide the types of logging and
    rollbacks that we desire. It is callable and works as a context handler to
    manage sessions when doing so.

    Attributes
    ----------
    logger : logging.Logger
        The logger object.
    maker : sqlalchemy.orm.session.Session
        A session maker.
    """

    def __init__(self, session_maker):
        """Create a new Session object.

        Parameter
        ---------
        session_maker : Session or sqlalchemy.orm.session.Session
            The session wrapper to use.
        """
        self.logger = logging.getLogger('core.Session')
        self.maker = session_maker
        if isinstance(session_maker, Session):
            self.maker = session_maker.maker

    @contextmanager
    def __call__(self):
        """Context handler for the session. This creates a new session and
        yields it. It will catch, log, and re raise any exceptions that occur.
        It will also commit and close all sessions. In the case of an error it
        will rollback the session.

        Yields
        ------
        session : Session
            A session that will be closed outside of the context.
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
