import logging
import collections as coll

from pymotifs.core.db import Session


class Base(object):
    """This is a simple utility class. Several things in core and outside of
    core need a basic class to inherit from that adds a logger, session
    handler, config, etc. This provides such a base class.
    """

    def __init__(self, config, session_maker, **kwargs):
        """Build a new Base object.

        :config: The config object to build with.
        :session_maker: The Session object to handle database connections.
        """

        self.config = coll.defaultdict(dict)
        self.config.update(config)
        self.session = Session(session_maker)
        self.name = self.__class__.__module__.replace('pymotifs.', '')
        self.logger = logging.getLogger(self.name)
