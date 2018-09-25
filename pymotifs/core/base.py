"""This module contains basic data structure that are needed by the entire
pipeline.
"""

import logging
import collections as coll

from pymotifs.core.db import Session


class Base(object):
    """This is a simple utility class. Several things in core and outside of
    core need a basic class to inherit from that adds a logger, session
    handler, config, etc. This provides such a base class.

    Attributes
    ----------
    config : dict
        The configuration object
    session : pymotifs.core.db.Session
        The session wrapper to use
    name : str
        The name of this stage. It will be the __name__ of the module but
        without the pymotifs component.
    logger : logging.Logger
        The logger to use.
    """

    def __init__(self, config, session_maker, **kwargs):
        """Build a new Base object.

        Parameters
        ----------
        config : dict
            The config object to build with.
        session_maker : pymotifs.core.db.Session
            The Session object to handle database connections.
        """

        self.config = coll.defaultdict(dict)
        self.config.update(config)
        self.session = Session(session_maker)
        self.name = self.__class__.__module__.replace('pymotifs.', '')
        self.logger = logging.getLogger(self.name)

    def _create(self, klass):
        return klass(self.config, self.session)
