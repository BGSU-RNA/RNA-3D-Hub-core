import os
import logging

logger = logging.getLogger(__name__)

try:
    from mlabwrap import mlab
except ImportError as err:
    logger.warning("No mlabwrap present, matlab will be skipped")
    logger.exception(err)
except Exception as err:
    logger.error("Could not import matlab")
    logger.exception(err)

from pymotifs.core.exceptions import Skip

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


class MatlabFailed(Exception):
    """An exception meant to be used if matlab commands fail for some reason.
    """
    pass


class Matlab(object):
    """A simple wrapper around mlab. This is useful because it sets the root as
    well as calling the setup function before running matlab.
    """

    def __init__(self, root):
        self.logger = logging.getLogger('core.Matlab')
        self.mlab = None
        self._root = root

    def __del__(self):
        if self.mlab:
            del self.mlab

    def __startup__(self):
        self.logger.debug('Starting up matlab')
        if 'mlab' not in globals():
            raise Skip("No matlab around, skipping")
        self.mlab = mlab
        os.chdir(BASE)
        # self.mlab._autosync_dirs = False
        self.mlab.setup()
        os.chdir(self._root)
        self.logger.debug('Matlab started')

    def __getattr__(self, key):
        if self.mlab is None:
            os.chdir(self._root)
            self.__startup__()
        self.logger.debug("Running %s", key)

        attr = getattr(self.mlab, key)

        def func(*args, **kwargs):
            corrected = []
            for arg in args:
                if isinstance(arg, basestring):
                    corrected.append(str(arg))
                else:
                    corrected.append(arg)
            return attr(*corrected, **kwargs)

        return func
