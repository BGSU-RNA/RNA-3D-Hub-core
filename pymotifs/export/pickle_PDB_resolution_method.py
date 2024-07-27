"""
Module for export of resolution and methods of pdb files
"""

from collections import defaultdict
import os
import pickle

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import DATA_FILE_DIRECTORY
from pymotifs.pdbs.info import Loader as InfoLoader


class Exporter(core.Loader):
    """

    """

    # General Setup
    compressed = False
    mark = False
    dependencies = set([InfoLoader])


    def has_data(self, *args, **kwargs):
        filename = self.filename()
        return False
        if os.path.exists(filename) is True:
            return True
        return False


    def remove():
        pass


    def filename(self):
        """
        Create the filename for the given PDB.

        Parameters
        ----------
        pdb :

        Returns
        -------
        filename : str
            The path to write to.
        """

        filename = 'PDB_resolution_method.pickle'

        self.logger.info("filename: filename: %s" % filename)

        return os.path.join(DATA_FILE_DIRECTORY,filename)


    def to_process(self, pdbs, **kwargs):
        """
        Ignore the pdbs input.
        The return value is just a list value.
        We do this because we only want to run this script one time even if we have a lot of pdbs.
        """

        if len(pdbs) < 100:
            raise core.Skip("Too few pdb files being processed to write PDB_resolution_method.pickle")

        return ['1']


    def data(self, pdb, **kwargs):
        """
        Look up all the existing pdbs to process.  Ignores the pdb input.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.pdb_id,
                            mod.ChainInfo.chain_name,
                            mod.ChainInfo.entity_macromolecule_type,
                            mod.PdbInfo.resolution,
                            mod.PdbInfo.experimental_technique).\
                            join(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id)
        result = defaultdict(dict)
        dna_long = ['Polydeoxyribonucleotide (DNA)','polydeoxyribonucleotide']
        rna_long = ['Polyribonucleotide (RNA)','polyribonucleotide']
        hybrid_long = ['polydeoxyribonucleotide/polyribonucleotide hybrid','DNA/RNA Hybrid']
        entity_type_check = dna_long + rna_long + hybrid_long + ['Polypeptide(L)','Peptide nucleic acid']

        for row in query:
            if not result.get(row.pdb_id):
                result[row.pdb_id]['chains'] = {}
            result[row.pdb_id]['resolution'] = row.resolution
            result[row.pdb_id]['method'] = row.experimental_technique

            if row.entity_macromolecule_type == 'Polypeptide(L)':
                if result[row.pdb_id]['chains'].get('protein'):
                    result[row.pdb_id]['chains']['protein'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['protein'] = []
                    result[row.pdb_id]['chains']['protein'].append(row.chain_name)
            elif (row.entity_macromolecule_type in dna_long):
                if result[row.pdb_id]['chains'].get('DNA'):
                    result[row.pdb_id]['chains']['DNA'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['DNA'] = []
                    result[row.pdb_id]['chains']['DNA'].append(row.chain_name)
            elif (row.entity_macromolecule_type in rna_long):
                if result[row.pdb_id]['chains'].get('RNA'):
                    result[row.pdb_id]['chains']['RNA'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['RNA'] = []
                    result[row.pdb_id]['chains']['RNA'].append(row.chain_name)
            elif (row.entity_macromolecule_type in hybrid_long):
                if result[row.pdb_id]['chains'].get('hybrid'):
                    result[row.pdb_id]['chains']['hybrid'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['hybrid'] = []
                    result[row.pdb_id]['chains']['hybrid'].append(row.chain_name)
            elif (row.entity_macromolecule_type == 'Peptide nucleic acid'):
                if result[row.pdb_id]['chains'].get('PNA'):
                    result[row.pdb_id]['chains']['PNA'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['PNA'] = []
                    result[row.pdb_id]['chains']['PNA'].append(row.chain_name)
            elif (row.entity_macromolecule_type not in entity_type_check):
                if result[row.pdb_id]['chains'].get(row.entity_macromolecule_type):
                    result[row.pdb_id]['chains'][row.entity_macromolecule_type].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains'][row.entity_macromolecule_type] = []
                    result[row.pdb_id]['chains'][row.entity_macromolecule_type].append(row.chain_name)

        return result


    def process(self, pdb, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        filename = self.filename()

        pinfo = self.data(pdb)

        with open(filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

