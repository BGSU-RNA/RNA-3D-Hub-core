import sys
import zlib
import logging
import traceback
import xml.etree.ElementTree as ET

from MotifAtlasBaseClass import MotifAtlasBaseClass

from models import NtQuality
from utils import FTPFetchHelper
from utils import DatabaseHelper

logger = logging.getLogger(__name__)


class FileHelper(object):
    path = 'pub/pdb/validation_reports/{short}/{pdb}/{pdb}_validation.xml.gz'

    def __call__(self, pdb):
        if not pdb:
            raise ValueError("Must give a pdb")

        return self.path.format(pdb=pdb, short=pdb[1:3])


class Parser(object):
    """A class to parse the results of getting the quality file.
    """

    def __call__(self, gz_content):
        content = zlib.decompress(gz_content)
        tree = ET.fromstring(content)
        root = tree.getroot()
        pass


class Loader(MotifAtlasBaseClass, DatabaseHelper):
    def __init__(self, session_builder):
        self.fetcher = FTPFetchHelper('ftp.wwpdb.org', parser=Parser())
        self.finder = FileHelper()
        MotifAtlasBaseClass.__init__(self)
        DatabaseHelper.__init__(self, session_builder)

    def data(self, pdb):
        filename = self.finder(pdb)
        logger.info("Using filename %s for pdb %s ", filename, pdb)
        for entry in self.fetcher(filename):
            yield NtQuality(**entry)

    def __call__(self, pdbs, **kwargs):
        for pdb in pdbs:
            logger.info("Fetching quality data for %s", pdb)
            try:
                data = self.data(pdb, **kwargs)
                self.store(data)
            except:
                logger.error("Error raised in getting quality for %s", pdb)
                logger.error(traceback.format_exc(sys.exc_info()))
