import inspect

from pymotifs.utils import known_subclasses

from .core import Representative

from .using_structure import *
from .using_quality import *


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
    for key, value in globals().items():
        if inspect.isclass(value) and issubclass(value, Representative) and \
                getattr(value, 'method', None) == name:
            return value
    raise ValueError("Unknown method %s" % name)
