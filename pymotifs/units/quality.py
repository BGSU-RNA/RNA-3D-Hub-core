import gzip
import hashlib
import cStringIO as sio
import collections as coll
import xml.etree.ElementTree as ET

try:
    from elementtree import SimpleXMLTreeBuilder
    ET.XMLTreeBuilder = SimpleXMLTreeBuilder.TreeBuilder
except:
    pass

import pymotifs.utils as ut
import pymotifs.core as core
from pymotifs import models as mod
# import UnitInfo as Unit
# from pymotifs.models import UnitQuality as Quality
from pymotifs.units.info import Loader as InfoLoader

from fr3d.unit_ids import encode


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
        self.generator = encode
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
        insertion_code = attributes['icode'].strip()
        if insertion_code == '':
            insertion_code = None

        alt_id = attributes.get('altcode', '').strip()
        if alt_id == '':
            alt_id = None

        return {
            'pdb': pdb,
            'model': int(attributes['model']),
            'chain': attributes['chain'],
            'component_number': int(attributes['resnum']),
            'component_id': attributes['resname'],
            'insertion_code': insertion_code,
            'alt_id': alt_id,
        }


class Loader(core.SimpleLoader):
    dependencies = set([InfoLoader])
    allow_no_data = True
    table = mod.UnitQuality

    def __init__(self, *args):
        super(Loader, self).__init__(*args)
        self.fetcher = ut.FTPFetchHelper('ftp.wwpdb.org')
        self.finder = FileHelper()

    def query(self, session, pdb):
        return session.query(mod.UnitQuality).\
            join(mod.UnitInfo, mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
            filter(mod.UnitInfo.pdb_id == pdb)

    def mapping(self, pdb):
        mapping = coll.defaultdict(list)
        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter_by(pdb_id=pdb)

            for result in query:
                key = (result.chain, result.number, result.ins_code,
                       result.alt_id)
                mapping[key].append(result.unit_id)

        return mapping

    def as_quality(self, entry, mapping):
        key = (entry['id']['chain'], entry['id']['component_number'],
               entry['id']['insertion_code'], entry['id']['alt_id'])

        for unit_id in mapping[key]:
            yield {
                'unit_id': unit_id,
                'real_space_r': entry.get('real_space_r'),
                'density_correlation': entry.get('density_correlation'),
                'z_score': entry.get('z_score')
            }

    def data(self, pdb, **kwargs):
        filename = self.finder(pdb)
        try:
            response = self.fetcher(filename)
        except ut.RetryFailedException:
            raise core.Skip("Could not download data for %s" % pdb)

        parser = Parser(response)
        if not parser.has_rsr() and not parser.has_dcc():
            raise core.Skip("No RsR found for %s" % pdb)

        mapping = self.mapping(pdb)
        data = []
        for entry in parser.nts():
            data.extend(quality for quality in self.as_quality(entry, mapping))
        return data
