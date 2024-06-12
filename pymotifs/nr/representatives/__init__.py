import inspect

from pymotifs.utils import known_subclasses

from .core import Representative

from .using_structure import *
from .using_quality import *
import time
import logging


def known():
    """
    List all the known methods to select representatives.
    """

    finders = []
    for subclass in known_subclasses(Representative, globals()):
        finders.append((subclass.method, subclass))
    return finders


def fetch(name):
    """
    Get the class that implements the given method name.

    Parameters
    ----------
    name : str
        Method name to use

    Returns
    -------
    finder : class
        A class that implements the given method.
    """
    logger = logging.getLogger(__name__)

    if not hasattr(fetch, 'call_count'):
        fetch.call_count = 0
    if not hasattr(fetch, 'total_time'):
        fetch.total_time = 0
    # Record the start time
    start_time = time.time()

    for key, value in globals().items():
        if inspect.isclass(value) and issubclass(value, Representative) and \
                getattr(value, 'method', None) == name:

            # Record the end time ... but it doesn't really work right
            # end_time = time.time()
            fetch.call_count += 1
            # fetch.total_time += (end_time - start_time)
            logger.info("Updated representative set number %d" % (fetch.call_count)
            # logger.info("show the name of the representative set: %s", name)
            # logger.info("show the representative info: %s", value)
            return value

    # # Record the end time
    # end_time = time.time()
    # fetch.call_count += 1
    # fetch.total_time += (end_time - start_time)
    # self.logger.info("Record the current running times and time for searching representatives. Running %s times, and cost %s",fetch.call_count,fetch.total_time)

    raise ValueError("Unknown method %s" % name)
