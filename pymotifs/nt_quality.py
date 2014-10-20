import gzip
import hashlib
import logging
import cStringIO as sio
import xml.etree.ElementTree as ET

import utils as ut
import core
from models import NtQuality
from models import PdbUnitIdCorrespondence as Unit

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
                data['real_space_r'] =  float(residue.attrib['rsr'])

            if 'DCC' in residue.attrib:
                pass

            if data:
                data['unit_id'] =  self._unit_id(pdb, residue.attrib)
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


class Loader(core.Loader):
    name = 'nt_quality'
    update_gap = False

    def __init__(self, config, session_builder):
        self.fetcher = ut.FTPFetchHelper('ftp.wwpdb.org')
        self.finder = FileHelper()
        super(Loader, self).__init__(self, config, session_builder)

    def remove(self, pdb):
        with self.session() as session:
            session.query(NtQuality).\
                join(Unit, Unit.unit_id == NtQuality.unit_id).\
                filter(Unit.pdb == pdb).\
                delete()

    def has_data(self, pdb):
        with self.session() as session:
            count = session.query(NtQuality).\
                join(Unit, Unit.unit_id == NtQuality.unit_id).\
                filter(Unit.pdb == pdb).\
                count()
            return bool(count)

    def data(self, pdb, **kwargs):
        filename = self.finder(pdb)
        logger.info("Using filename %s for pdb %s ", filename, pdb)
        response = self.fetcher(filename)
        parser = Parser(response)
        if not parser.has_rsr() and not parser.has_dcc():
            logger.info("No RsR found for %s", pdb)
            return

        for entry in parser.nts():
            yield NtQuality(**entry)

if __name__ == '__main__':
    ut.main(Loader)
