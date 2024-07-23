"""A module containing general utility functions. This are common utilities
that show up across the pipeline.
"""

import os
import gzip
import inspect
import logging
from ftplib import FTP
import itertools as it
import cStringIO as sio
import collections as coll

import requests


"""Generic logger for all utilities."""
logger = logging.getLogger(__name__)


class MissingFileException(Exception):
    """This is raised when we can't find a file. For example a cif file for a
    pdb does not exist.
    """
    pass


class WebRequestFailed(Exception):
    """Raised when we have failed all attempts at getting a url.
    """
    pass


class EmptyResponse(Exception):
    """Raised when processing an empty response.
    """
    pass


class RetryFailedException(Exception):
    """Raised when all attempts at retrying something have failed.
    """
    pass


class GetAllRnaPdbsError(Exception):
    """Raised with we cannot get all RNA containing PDBS
    """
    pass


class GetCustomReportError(Exception):
    """Raised when we cannot get a custom report
    """
    pass


def grouper(n, iterable):
    """Group an iterable in chunks of a max size.

    Parameters
    ----------
    n : int
        The max size.
    iterable : iterable
        The iterable to group.

    Returns
    -------
    chunks : iterable
        A iterable with at most `n` entries per chunk.
    """

    iterator = iter(iterable)
    while True:
        chunk = tuple(it.islice(iterator, n))
        if not chunk:
            return
        yield chunk


def list_or_tuple(obj):
    """Detect if something is a list or a tuple. This is useful when flattening
    lists with flatten.

    Parameters
    ----------
    obj : object
        An object

    Returns
    -------
    result : bool
        True if the object is a list or tuple.
    """

    return isinstance(obj, (list, tuple))


def flatten(sequence, predicate=list_or_tuple):
    """Flatten deeply nested lists. This produces a generator which will
    flatten deepely nested lists. Modified from:
    https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch04s07.html

    Parameters
    ----------
    sequence: list or tuple or generator
        A list, tuple or generator to flatten.

    predicate : function, optional
        The predicate to use to detect if we should flatten the contents of
        each object in the sequence.

    Yields
    ------
    component : object
        The individual components from the nested iterator.
    """

    for obj in sequence:
        if predicate(obj):
            for entry in flatten(obj):
                yield entry
        else:
            yield obj


def row2dict(row):
    """Convert an sqlalchemy object into a dictionary. This should either be
    the result of querying or an object that is to be saved. It must have
    either a '__table__' property or a '_feilds' property.

    Parameters
    ----------
    row
        Object to convert.

    Returns
    -------
    result : dict
        A dictionary built from the contents of the object.
    """

    if row is None:
        return None
    d = {}
    if hasattr(row, '__table__'):
        for column in row.__table__.columns:
            value = getattr(row, column.name)
            d[column.name] = value
    elif hasattr(row, '_fields'):
        for column in row._fields:
            d[column] = getattr(row, column)
    return d


def result2dict(result):
    """Turn a result from a query into a dictionary

    Parameters
    ----------
    result : object
        A result from a sqlalchemy query.

    Returns
    -------
    data : dict
        A dictionary with the keys from the query.
    """

    d = {}
    for name in result.keys():
        d[name] = getattr(result, name)
    #     print("name:",name)
    #     print("result:",result)
    # print("result.keys():",result.keys())
    # print("d:",d)
    return d


def known_subclasses(base, glob):
    """Get the list of known subclasses from the given dictonary.
    """
    subclasses = []
    for key, value in glob.items():
        if inspect.isclass(value) and \
                issubclass(value, base) and \
                value is not base:
            subclasses.append(value)
    return subclasses


def known(config, pdb=False, cif=True, pdb1=False, cifatoms=True):
    """Find all known mmCIF, pdb and/or PDB1 files. Giving all flags means
    getting all PDB, mmCIF or PDB1 files. This will check if we have a copy of
    the file downloaded with the correct extension.

    Parameters
    ----------
    config : defaultdict
        The configuration object that species the 'fr3d_root'.
    pdb : bool, optional
        A flag if we should get known PDB files.
    cif : bool, optional
        A flag if we should get known mmCIF files.
    pdb1 : bool, optional
        A flag if we should get known pdb1 files.

    Yields
    ------
    filename : str
        A PDB id that is known.
    """

    path = os.path.join(config['locations']['fr3d_root'], 'PDBFiles')
    names = coll.defaultdict(dict)

    for filename in os.listdir(path):
        if not os.path.isfile(os.path.join(path, filename)):
            continue
        if 'exemplar' in filename.lower():
            continue
        name, ext = os.path.splitext(filename)
        ext = ext.replace('.', '')
        names[name][ext] = True

    for name, exts in names.items():
        if not pdb or (pdb and exts.get('pdb')):
            if not cif or (cif and exts.get('cif')):
                if not pdb1 or (pdb1 and exts.get('pdb1')):
                    if not cifatoms or (cifatoms and exts.get('cifatoms')):
                        yield name


class RetryHelper(object):
    """A base class for retrying an action. Some actions, like making an HTTP
    request can fail because of intermittent issues, like network failure.
    For things like this it makes sense to

    Attributes
    ----------
    retries : int
        The number of retries.
    allow_fail : bool
    log_last_only : bool
        Only log the last failure, not all the initial ones.
    """

    def __init__(self, retries=3, allow_fail=False, log_last_only=True):
        """Create a new Retry helper.
        """
        self.retries = retries
        self.allow_fail = allow_fail
        self.log_last_only = log_last_only

    def action(self, *args, **kwargs):
        """The action method that all retry helpers should implement.
        """
        pass

    def __call__(self, *args, **kwargs):
        """Attempt the action with the given number of retries.
        """

        for index in xrange(self.retries):
            try:
                return self.action(*args, **kwargs)
            except Exception as err:
                if not self.log_last_only or \
                        (self.log_last_only and index == (self.retries - 1)):

                    if not self.allow_fail:
                        logger.warning("Failed retry attempt #%s",
                                       str(index + 1))
                        logger.exception(err)

        if not self.allow_fail:
            logger.error("All attempts at retrying failed")
            raise RetryFailedException()


class WebRequestHelper(RetryHelper):
    """
    A class to help with making web requests. This deals with the retrying
    and making sure the response body is not empty. If the max number of
    retries is reached we raise an exception. In addition, all steps are
    logged.

    Attributes
    ----------
    allow_empty : bool
        Flag to indicate if we should allow an empty response.
    method : str
        HTTP method to use.
    parser : function
        A function to parse the response body with.
    """

    def __init__(self, allow_empty=False, method='get', parser=None, **kwargs):
        """Create a new WebRequest helper.

        Parameters
        ----------
        allow_empty : bool
            Flag to indicate if we should allow an empty response.
        method : str
            HTTP method to use.
        parser : function
            A function to parse the response body with.
        """

        self.method = method
        self.allow_empty = allow_empty
        self.parser = parser
        super(WebRequestHelper, self).__init__(**kwargs)

    def action(self, url, **kwargs):
        method = getattr(requests, self.method)
        logger.info("Sending request to %s", url)
        response = method(url, **kwargs)
        response.raise_for_status()

        if not self.allow_empty and not response.text and not self.allow_fail:
            logger.warning("Response body was empty. Retrying.")
            raise WebRequestFailed("Got empty response")
        if self.parser:
            return self.parser(response)
        return response.text


class StructureFileFinder(object):
    extension = None

    def __init__(self, config, extension=None, strict=True):
        self.config = config
        self.strict = strict
        self.location = os.path.join(self.config['locations']['fr3d_root'],
                                     'PDBFiles')
        if extension:
            self.extension = extension

        if not self.extension:
            raise ValueError("Must define an extension")

    def __call__(self, pdb):
        filename = os.path.join(self.location, pdb + '.' + self.extension)
        filename = os.path.realpath(os.path.abspath(filename))

        if not os.path.exists(filename):
            if not self.strict:
                return None

            msg = "Could not find CIF file for %s.  Expected at: %s"
            raise MissingFileException(msg % (pdb, filename))

        return filename


class CifFileFinder(StructureFileFinder):
    """This is a class that can find cif files. It looks into the configure
    PDBFiles directory for cif files. It will return the full path to the file,
    or raise an exception if it could not be found in strict mode. If not
    strict it returns None.
    """

    """File extension to find."""
    extension = 'cif'


class PickleFileFinder(StructureFileFinder):
    extension = 'pickle'


class FTPFetchHelper(RetryHelper):
    def __init__(self, uri, allow_empty=False, parser=None, **kwargs):
        self.parser = parser
        self.allow_empty = allow_empty
        logger.debug("Connecting to %s", uri)
        self.ftp = FTP(uri)
        self.ftp.login()
        logger.debug("Connected and logged in")
        super(FTPFetchHelper, self).__init__(**kwargs)

    def action(self, filename, **kwargs):
        logger.info("Attempting to get %s", filename)
        out = sio.StringIO()
        self.ftp.retrbinary("RETR %s" % filename, out.write)
        logger.info("Fetched %s", filename)
        text = out.getvalue()

        if not self.allow_empty and not text:
            logger.warning("Retrived content was empty.")
            raise EmptyResponse("Got empty content.")

        if self.parser:
            return self.parser(text)
        return text


class GzipFetchHelper(WebRequestHelper):
    def __init__(self, **kwargs):
        super(GzipFetchHelper, self).__init__(parser=self.parser)

    def parser(self, response):
        # logger.info('response is %s', response)
        # logger.info('response content is %s', response.content)  # prints binary data from .cif.gz!
        # logger.info('response is %s', response)

        fileobj = sio.StringIO(response.content)
        unzipped = gzip.GzipFile(fileobj=fileobj)
        return unzipped.read()
