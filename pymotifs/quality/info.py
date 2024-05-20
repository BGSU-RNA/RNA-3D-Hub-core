"""This will download information about
"""

import os
import gzip
import hashlib
import operator as op
import itertools as it
import functools as ft
import collections as coll
import xml.etree.ElementTree as ET

import sys

# safe for python 2 and 3
if sys.version_info[0] == 2:
    import cStringIO as sio
else:
    import io as sio

try:
    from elementtree import SimpleXMLTreeBuilder
    ET.XMLTreeBuilder = SimpleXMLTreeBuilder.TreeBuilder
except:
    pass

import pymotifs.utils as ut
import pymotifs.core as core
from pymotifs import models as mod

from pymotifs.units.info import Loader as InfoLoader
from pymotifs.quality.download import Loader as Downloader

from fr3d.unit_ids import encode


"""A `operator.itemgetter` for defining a key for mapping between unit ids and
entries in the validation report.
"""
as_key = op.itemgetter('chain', 'number', 'ins_code', 'alt_id')


class FileHelper(object):
    """Class to help computing the path to the validation report on PDB's FTP
    site.
    """

    """The template for paths to the validation report"""
    path = 'pub/pdb/validation_reports/{short}/{pdb}/{pdb}_validation.xml.gz'

    def __call__(self, pdb):
        """Compute the path to the validation report for the given file.

        Parameters
        ----------
        pdb : str
            The pdb id to get a report for.

        Returns
        -------
        path : str
            The path on the FTP server to use.
        """

        if not pdb:
            raise ValueError("Must give a pdb")

        return self.path.format(pdb=pdb.lower(), short=pdb[1:3].lower())


class Parser(object):
    """A class to parse the results of getting the quality file. Right now it
    only processes the RsR related data.

    Attributes
    ----------
    generator : fr3d.unit_ids.encode
        A function to generate unit ids from dictionaries.
    digest : str
        The md5 hash of the content
    root : ElementTree
        The parsed XML tree for processing.
    """

    """List of keys from the validation report to extract"""
    entry_keys = [
        'absolute-percentile-percent-RSRZ-outliers',
        'high-resol-relative-percentile-percent-RSRZ-outliers',
        'low-resol-relative-percentile-percent-RSRZ-outliers',
        'numPDBids-absolute-percentile-percent-RSRZ-outliers',
        'numPDBids-relative-percentile-percent-RSRZ-outliers',
    ]

    def __init__(self, gz_content):
        """Create a new `Parser` to parse the given gz_content. This parser
        will extract the residue level entries from the content.

        Parameters
        ----------
        gz_content : str
            A gzip'ed string of the file to parse.
        """

        filehandle = sio.StringIO(gz_content)
        content = gzip.GzipFile(fileobj=filehandle).read()
        md5 = hashlib.md5()
        md5.update(content)
        self.generator = encode
        self.digest = md5.hexdigest()
        self.root = ET.fromstring(content)

    def has_dcc(self):
        """Check if this report has DCC data.

        Returns
        -------
        has_dcc : bool
            True if this report has DCC data.
        """

        entry = self.root.find("Entry")
        return 'DCC_R' in entry.attrib

    def has_rsr(self):
        """Check if this report has RSR data.

        Returns
        -------
        True if this report has RSR data
        """

        entry = self.root.find("Entry")
        return 'absolute-percentile-percent-RSRZ-outliers' in entry.attrib

    def nts(self):
        """Get all nucleotide data from the parsed tree. This will extract all
        residue level quality data and produce an iterator over the resulting
        dictionaries. This will extract data for all entries, RNA, DNA, and all
        ligands.

        Yields
        ------
        nt : dict
            A dictionary of nt level data.
        """

        pdb = self.root.find("Entry").attrib['pdbid'].upper()
        for residue in self.root.findall("ModelledSubgroup"):
            data = {}
            if 'rsr' in residue.attrib:
                data['real_space_r'] = float(residue.attrib['rsr'])

            if 'rsrz' in residue.attrib:
                data['real_space_r_z_score'] = float(residue.attrib['rsrz'])

            if 'DCC' in residue.attrib:
                data['density_correlation'] = float(residue['DCC'])

            if data:
                data['id'] = self._unit_id(pdb, residue.attrib)
                yield data

    def _unit_id(self, pdb, attributes):
        """Compute a dictionary that can be used to find and compute unit ids.
        This does not create a real unit id string but instead creates a
        dictionary that can later be used to generate one. This does not add
        data like the symmetry operator as this cannot be determined from the
        validation report alone. That has to be inferred by finding all
        symmetry operators for units that this unit id data can refer to.

        Parameters
        ----------
        pdb : str
            The PDB Id
        attributes : dict
            Dictionary that includes 'model', 'chain', 'resnum', 'resname', and
            optionally, 'icode' and 'altcode'.

        Returns
        -------
        unit_id : dict
            A dictionary with keys, 'pdb', 'model', 'chain', component_number',
            'component_id', 'insertion_code', and 'alt_id'.
        """

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
            'number': int(attributes['resnum']),
            'component_id': attributes['resname'],
            'ins_code': insertion_code,
            'alt_id': alt_id,
        }


class Loader(core.SimpleLoader):
    """The loader to fetch and store quality data for structures.
    """

    dependencies = set([InfoLoader, Downloader])

    def to_process(self, pdbs, **kwargs):
        empty = set(f for f in self.known() if os.stat(f).st_size == 0)
        return sorted(set(pdbs) - empty)

    def query(self, session, pdb):
        """Generate a query to find all entries in units_quality for the given
        PDB id.

        Attributes
        ----------
        session : Session
            The `Session` to use.

        pdb : str
            The PDB id to use.

        Returns
        -------
        query : Query
            Returns an SqlAlchemy query for all entires in units_quality with
            the given PDB id.
        """

        return session.query(mod.UnitQuality).\
            join(mod.UnitInfo,
                 mod.UnitInfo.unit_id == mod.UnitQuality.unit_id).\
            filter(mod.UnitInfo.pdb_id == pdb)

    def mapping(self, pdb):
        """Create a dictionary that maps from data produced by `as_key` to unit
        ids that are in the database. This will lookup all unit ids in the
        database and create the required mapping.

        Parameters
        ----------
        pdb : str
            The pdb id to look up a mapping for

        Returns
        -------
        mapping : dict
            The mapping dictionary to use.
        """

        mapping = coll.defaultdict(list)
        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter_by(pdb_id=pdb)

            for result in query:
                key = as_key(ut.row2dict(result))
                mapping[key].append(result.unit_id)

        return mapping

    def as_quality(self, mapping, entry):
        """Convert an entry from the parser into a form suitable for writing to
        the units_quality table. Since some entries from the parser expand to
        more than one unit due to symmetry operators this will produce an
        iterator that may have more than 1 value.

        Parameters
        ----------
        mapping : dict
            The mapping as produced by `mapping`.
        entry : dict
            A dictionary from `Parser.nts`.

        Yields
        ------
        entry : dict
            A dictionary of 'unit_id', 'real_space_r', 'density_correlation',
            'real_space_r_z_score'.
        """

        key = as_key(entry['id'])

        if not mapping[key]:
            raise core.InvalidState("Could not find unit id for %s" % entry)

        for unit_id in mapping[key]:
            yield {
                'unit_id': unit_id,
                'real_space_r': entry.get('real_space_r'),
                'density_correlation': entry.get('density_correlation'),
                'real_space_r_z_score': entry.get('real_space_r_z_score')
            }

    def data(self, pdb, **kwargs):
        """Compute the quality assignments for residues in the structure. This
        will fetch the validation report from PDB and convert the entries there
        into forms suitable to write to the database. If the report has no RSR
        or DCC data then a `core.Skip` exception will be raised.

        Parameters
        ----------
        pdb : str
            The pdb id to use.

        Returns
        -------
        data : iterable
            An iterable of a quality assignments to store in the database.
        """

        with open(self.filename(pdb), 'rb') as raw:
            parser = Parser(raw.read())

        if not parser.has_rsr() and not parser.has_dcc():
            raise core.Skip("No RsR found for %s" % pdb)

        mapping = self.mapping(pdb)
        as_quality = ft.partial(self.as_quality, mapping)
        data = it.imap(as_quality, parser.nts())
        data = it.chain.from_iterable(data)
        return it.imap(lambda d: mod.UnitQuality(**d), data)

