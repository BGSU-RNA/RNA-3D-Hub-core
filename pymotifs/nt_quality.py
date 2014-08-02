import sys
import gzip
import hashlib
import logging
import traceback
import cStringIO as sio
import xml.etree.ElementTree as ET

from MotifAtlasBaseClass import MotifAtlasBaseClass

import utils as ut
from models import NtQuality

ut.append_libs()
from rnastructure.util import unit_ids as uids


logger = logging.getLogger(__name__)


class FileHelper(object):
    path = 'pub/pdb/validation_reports/{short}/{pdb}/{pdb}_validation.xml.gz'

    def __call__(self, pdb):
        if not pdb:
            raise ValueError("Must give a pdb")

        return self.path.format(pdb=pdb.lower(), short=pdb[1:3].lower())


class Parser(object):
    """A class to parse the results of getting the quality file. Right now it
    only processes the RsR related data.
    """
    entry_keys = [
        'absolute-percentile-percent-RSRZ-outliers',
        'high-resol-relative-percentile-percent-RSRZ-outliers',
        'low-resol-relative-percentile-percent-RSRZ-outliers',
        'numPDBids-absolute-percentile-percent-RSRZ-outliers',
        'numPDBids-relative-percentile-percent-RSRZ-outliers',
    ]

    def __init__(self, gz_content):
        filehandle = sio.StringIO(gz_content)
        content = gzip.GzipFile(fileobj=filehandle).read()
        md5 = hashlib.md5()
        md5.update(content)
        self.generator = uids.UnitIdGenerator()
        self.digest = md5.hexdigest()
        self.root = ET.fromstring(content)
        print(self.root)  # .getroot()

    def has_rsr(self):
        entry = self.root.find("Entry")
        return 'absolute-percentile-percent-RSRZ-outliers' in entry.attrib

    def nts(self):
        pdb = self.root.find("Entry").attrib['pdbid'].upper()
        for residue in self.root.findall("ModelledSubgroup"):
            if 'rsr' in residue.attrib and 'rsrz' in residue.attrib:
                yield {
                    'id': self._unit_id(pdb, residue.attrib),
                    'rsr': float(residue.attrib['rsr']),
                    'rsrz': float(residue.attrib['rsrz'])
                }

    def _unit_id(self, pdb, attributes):
        return self.generator({
            'pdb': pdb,
            'model': int(attributes['model']),
            'chain': attributes['chain'],
            'number': int(attributes['resnum']),
            'residue': attributes['resname'],
            'insertion_code': attributes['icode'],
            'symmetry_operator': '*'
        })


class Loader(MotifAtlasBaseClass, ut.DatabaseHelper):
    def __init__(self, session_builder):
        self.fetcher = ut.FTPFetchHelper('ftp.wwpdb.org')
        self.finder = FileHelper()
        MotifAtlasBaseClass.__init__(self)
        ut.DatabaseHelper.__init__(self, session_builder)

    def data(self, pdb):
        filename = self.finder(pdb)
        logger.info("Using filename %s for pdb %s ", filename, pdb)
        response = self.fetcher(filename)
        parser = Parser(response)
        if not parser.has_rsr():
            logging.info("No RsR found for %s", pdb)
            return

        for entry in parser.nts():
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
