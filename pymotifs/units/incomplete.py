"""This is a loader to store data on which residues in a structure are
incomplete. Units can be incomplete in that they have atoms (or entire
residues) that have 0 occpancy. This is an issue that shows up in grouping
loops later.
"""

import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod

from pymotifs.units.info import Loader as InfoLoader


class Entry(coll.namedtuple('Entry', ['pdb_id', 'model', 'chain', 'number',
                                      'unit', 'alt_id', 'ins_code'])):
    """A useful tuple for creating mappings and handling incomplete/missing
    units. It does not include symmetry operators to allow for mapping from
    these entries to all unit ids which will include symmetry operators.

    Attributes
    ----------
    model : int
        The model number
    chain : str
        The chain name
    number : int
        The unit number
    sequence : str
        The sequence of the unit
    alt_id : str or None
        The alt id
    insertion : str or None
        The insertation code
    """
    pass


class Loader(core.SimpleLoader):
    merge = True
    allow_no_data = True
    dependencies = set([InfoLoader])

    def query(self, session, pdb):
        """Build a query to

        Parameters
        ----------
        session : pymotifs.core.Session
            The session to use.

        pdb : str
            The pdb to query for.
        """
        return session.query(mod.UnitIncomplete).\
            filter(mod.UnitIncomplete.pdb_id == pdb)

    def missing_keys(self, pdb):
        """Determine the unit ids in a file that are missing/unobserved. This
        parses the mmCIF file for the given PDB and examines the
        'pdbx_unobs_or_zero_occ_residues' and 'pdbx_unobs_or_zero_occ_atoms'
        data blocks for this information.

        Parameters
        ----------
        pdb : str
            The PDB id.

        Returns
        -------
        missing : set
            A set of Entry tuples
        """
        cif = self.cif(pdb)
        total = it.chain(getattr(cif, 'pdbx_unobs_or_zero_occ_residues', []),
                         getattr(cif, 'pdbx_unobs_or_zero_occ_atoms', []))
        missing = set()
        for row in total:
            model = int(row['PDB_model_num'])
            chain = row['auth_asym_id']
            num = int(row['auth_seq_id'])
            seq = row['auth_comp_id']
            ins = None
            if row['PDB_ins_code'] != '?':
                ins = row['PDB_ins_code']

            alt_id = None
            if 'label_alt_id' in row and row['label_alt_id'] != '?':
                alt_id = row['label_alt_id']

            key = Entry(pdb, model, chain, num, seq, alt_id, ins)
            missing.add(key)
        return missing

    def data(self, pdb, **kwargs):
        """Determine the unit ids from a given structure that are incomplete
        and create data for writing to the database. This will only create
        entries for things that are incomplete, not missing because of database
        constraints. The unit_incomplete table requries that the unit_id
        exist in unit_info as well, and missing ones will not (by definition).


        Parameters
        ----------
        pdb : str
            The PDB id.

        Returns
        -------
        data : list
            A list of pymotifs.model.UnitIncomplete entries for all incomplete
            unit ids.
        """
        data = []
        for key in self.missing_keys(pdb):
            data.append(mod.UnitIncomplete(**key._asdict()))
        return data
