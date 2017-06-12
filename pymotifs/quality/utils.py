"""This module contains generic utiltiy functions and classes for working with
validation reports.
"""

import os
import copy
import gzip
import hashlib
import operator as op
import cStringIO as sio
import collections as coll
import xml.etree.ElementTree as ET

from pymotifs import core
from pymotifs import utils as ut
from pymotifs import models as mod
from pymotifs.utils import renaming as rn

from fr3d.unit_ids import encode


def as_key(entry):
    """A key to use for identify both unit ids and entries in the quality
    reports. The key should be considered opaque, but it will add a model entry
    even if the given data does not have one. If it does not the model will
    default to 1.

    :param dict entry: The entry to compute a key for.
    :returns: A tuple, that is a unique key for the entry.
    """

    simple = op.itemgetter('chain', 'number', 'ins_code', 'alt_id')
    current = list(simple(entry))
    current.append(entry.get('model', 1))
    return tuple(current)


class Utils(core.Base):
    """A set of a utilities for dealing with quality data.
    """

    def filename(self, pdb, compressed=True):
        """Compute the filename for the given PDB id.

        Parameters
        ----------
        pdb : str
            The PDB id to use.
        compressed : bool, True
            Flag if we should create a filename for the compressed report.

        Returns
        -------
        filename : str
            Filename for the given PDB.
        """
        base = self.config['locations']['quality_reports']
        ext = '.xml'
        if compressed:
            ext += '.gz'
        return os.path.join(base, pdb + ext)

    def known(self, has_data=None):
        """Get a list of all known PDBs that have downloaded data. By default this
        will list all PDBs, even those that have no quality data. However, this is
        changable using the has_data flag. Setting it to True will filter the list
        to only those that have data, while False, will be for those that were
        downloaded but have no data.

        Parameters
        ----------
        config : dict
            The configuration dictonary.

        has_data : bool, None
            A flag to indicate if we should require that files be non-empty (True),
            empty (False), or either (None).
        """
        pdbs = []
        dirname = self.config['locations']['quality_reports']
        for basename in os.listdir(dirname):
            filename = os.path.join(dirname, basename)
            if has_data is not None:
                if bool(os.stat(filename).st_size) != has_data:
                    continue
            pdbs.append(basename[0:4])
        return pdbs

    def has_no_data(self, pdb):
        """Check if the validation report for the given pdb is empty or not. If it
        is empty then this means there was no validation report. However, if the
        file was never attempted to be downloaded, this will fail.

        Parameters
        ----------
        config : dict
            The configuration dictonary

        pdb : str
            The PDB to use

        Returns
        -------
        has_no_data : bool
            True if the file for the PDB is empty and thus has no data.
        """
        name = self.filename(pdb)
        if not os.path.exists(name):
            return True
        return os.stat(name).st_size == 0

    def unit_mapping(self, pdb):
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

    unit_renamer = rn.Renamer(
        real_space_r=rn.rename('rsr', rn.maybe_float),
        real_space_r_z_score=rn.rename('rsrz', rn.maybe_float),
        density_correlation=rn.rename('DCC', rn.maybe_float),
    )

    unit_id_renamer = rn.Renamer(
        rn.transform('model', int),
        rn.transform('chain', str),
        component_id=rn.rename('resname', str),
        number=rn.rename('resnum', int),
        ins_code=rn.rename('icode', rn.maybe_str, strip=True),
        alt_id=rn.rename('altcode', rn.maybe_str, strip=True),
    )

    structure_renamer = rn.Renamer(
        rn.with_dashes('percent-RSRZ-outliers', rn.maybe_float),
        rn.with_dashes('absolute-percentile-percent-RSRZ-outliers', rn.maybe_float),
        rn.with_dashes('relative-percentile-percent-RSRZ-outliers', rn.maybe_float),
        rn.with_dashes('clashscore', rn.maybe_float),
        rn.with_dashes('relative-percentile-clashscore', rn.maybe_float),
        rn.with_dashes('absolute-percentile-clashscore', rn.maybe_float),
        rn.with_dashes('percent-rota-outliers', rn.maybe_float),
        rn.with_dashes('absolute-percentile-percent-rota-outliers', rn.maybe_float),
        rn.with_dashes('relative-percentile-percent-rota-outliers', rn.maybe_float),
        pdb_id=rn.rename('pdbid', rn.maybe_str),
    )

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

    def entity(self):
        """Get the entity level anotations.

        Returns
        -------
        entity : dict
            A dictonary of mappings for all attributes on the entity entry.
            The keys and values will all be strings.
        """
        data = self.structure_renamer(self.root.find("Entry").attrib)
        data['md5'] = self.digest
        return data

    def has_dcc(self):
        """Check if this report has DCC data.

        Returns
        -------
        has_dcc : bool
            True if this report has DCC data.
        """
        return 'DCC_R' in self.root.find("Entry").attrib

    def has_rsr(self):
        """Check if this report has RSR data.

        Returns
        -------
        True if this report has RSR data
        """
        entry = self.root.find("Entry")
        return 'absolute-percentile-percent-RSRZ-outliers' in entry.attrib

    def nts(self, mapping):
        """Get all nucleotide data from the parsed tree. This will extract all
        residue level quality data and produce an iterator over the resulting
        dictionaries. This will extract data for all entries, RNA, DNA, and all
        ligands.

        Parameters
        ----------
        mapping : dict
            A dictonary as from `Utils.mapping`. This will be used to map to
            known unit ids from entries in the validation file.

        Yields
        ------
        nt : dict
            A dictionary of nt level data.
        """
        for residue in self.root.findall("ModelledSubgroup"):
            data = self.unit_renamer(residue.attrib, skip_missing=True)
            if not data:
                continue

            uid = as_key(self.unit_id_renamer(residue.attrib))
            if uid not in mapping:
                raise core.InvalidState("Could not find unit id for %s" % data)

            if not mapping[uid]:
                raise core.InvalidState("No unit ids known for %s", uid)

            for pk in mapping[uid]:
                d = copy.deepcopy(data)
                d['id'] = pk
                yield d
