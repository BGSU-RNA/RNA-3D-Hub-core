import gzip
import hashlib
import cStringIO as sio
import xml.etree.ElementTree as ET

import pymotifs.utils as ut
import pymotifs.core as core
from pymotifs.models import UnitInfo as Unit
from pymotifs.models import UnitQuality as Quality

from rnastructure.util import unit_ids as uids


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

    def has_dcc(self):
        entry = self.root.find("Entry")
        return 'DCC_R' in entry.attrib

    def has_rsr(self):
        entry = self.root.find("Entry")
        return 'absolute-percentile-percent-RSRZ-outliers' in entry.attrib

    def nts(self):
        pdb = self.root.find("Entry").attrib['pdbid'].upper()
        for residue in self.root.findall("ModelledSubgroup"):
            data = {}
            if 'rsr' in residue.attrib:
                data['real_space_r'] = float(residue.attrib['rsr'])

            if 'DCC' in residue.attrib:
                pass

            if data:
                data['id'] = self._unit_id(pdb, residue.attrib)
                yield data

    def _unit_id(self, pdb, attributes):
        return self.generator({
            'pdb': pdb,
            'model': int(attributes['model']),
            'chain': attributes['chain'],
            'number': int(attributes['resnum']),
            'residue': attributes['resname'],
            'insertion_code': attributes['icode'],
        })


class Loader(core.SimpleLoader):

    def __init__(self, *args):
        super(Loader, self).__init__(*args)
        self.fetcher = ut.FTPFetchHelper('ftp.wwpdb.org')
        self.finder = FileHelper()

    def query(self, session, pdb):
        return session.query(Quality).\
            join(Unit, Unit.id == Quality.id).\
            filter(Unit.pdb_id == pdb)

    def data(self, pdb, **kwargs):
        filename = self.finder(pdb)
        response = self.fetcher(filename)
        parser = Parser(response)
        if not parser.has_rsr() and not parser.has_dcc():
            self.logger.info("No RsR found for %s", pdb)
            return

        return [Quality(**entry) for entry in parser.nts()]
