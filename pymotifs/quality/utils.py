import os
import gzip
import hashlib
import cStringIO as sio
import xml.etree.ElementTree as ET

from fr3d.unit_ids import encode


def filename(config, pdb, compressed=True):
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
    base = config['locations']['quality_reports']
    ext = '.xml'
    if compressed:
        ext += '.gz'
    return os.path.join(base, pdb + ext)


def known(config, has_data=None):
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
    for filename in os.path.listdir(config):
        if has_data is not None:
            if os.stat(filename).st_size != has_data:
                continue
        pdbs.append(os.path.basename(filename))
    return pdbs


def has_no_data(config, pdb):
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
    name = filename(config, pdb)
    return os.stat(name).st_size == 0


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

    def entity(self):
        """Get the entity level anotations.

        Returns
        -------
        entity : dict
            A dictonary of mappings for all attributes on the entity entry.
            The keys and values will all be strings.
        """
        return self.root.find("Entry").attrib

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
                data['z_score'] = float(residue.attrib['rsrz'])

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
