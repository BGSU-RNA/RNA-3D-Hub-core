"""
Module for export of resolution, chains, and experimental methods of existing pdb files
"""

import numpy as np
import os
import pickle

from pymotifs import core
from pymotifs import models as mod

from pymotifs.pdbs.info import Loader as InfoLoader

from collections import defaultdict

from os import path

from sqlalchemy import and_
from sqlalchemy.orm import aliased
from sqlalchemy.sql import select
from sqlalchemy.sql import union


class Exporter(core.Loader):
    """Export pairs data in pickle format.
    """

    # General Setup
    compressed = False
    mark = False
    dependencies = set([InfoLoader])
    #dependencies = set()

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

        # filename = 'PDB_resolution_method.pickle'
        filename = 'NA_datafile.pickle'

        self.logger.info("filename: %s" % filename)

        return os.path.join("/var/www/html/units",filename)

    def to_process(self, pdbs, **kwargs):
        """
        Ignore the pdbs input.
        The return value is just a list value.
        We do this because we only want to run this script one time even if we have a lot of pdbs.
        """

        if len(pdbs) < 100:
            raise core.Skip("Too few pdb files being processed to write NA_datafile.pickle")

        return ['1']


    def model_query(self, **kwargs):
        """
        Map PDB id to distinct model numbers
        """
        with self.session() as session:
            query = session.query(mod.UnitInfo.pdb_id,mod.UnitInfo.model).distinct()
        result = defaultdict(set)
        for row in query:
            try:
                result[row.pdb_id].add(int(str(row.model).replace('L','')))
            except:
                result[row.pdb_id].add(str(row.model).replace('L',''))
                if not row.pdb_id == 'XXXX':
                    self.logger.info('The database does not have the model type for %s' % row.pdb_id)
        return result


    def sym_op(self, **kwargs):
        """
        Map PDB id and chain to distinct symmetry operators
        """
        with self.session() as session:
            query = session.query(mod.UnitInfo.pdb_id,mod.UnitInfo.chain,mod.UnitInfo.sym_op).distinct()
        result = defaultdict(set)
        for row in query:
            result[str(row.pdb_id)+"_"+str(row.chain)].add(row.sym_op)
        return result

    def model_sym_op_query(self, **kwargs):
        """
        Map PDB id and chain to models and symmetry operators
        """

        pdb_id_to_models = defaultdict(set)
        pdb_chain_to_sym_op = defaultdict(set)

        with self.session() as session:
            query = session.query(mod.UnitInfo.pdb_id,mod.UnitInfo.model,mod.UnitInfo.chain,mod.UnitInfo.sym_op).distinct()

            for row in query:
                try:
                    model = int(str(row.model).replace('L',''))
                    pdb_id_to_models[row.pdb_id].add(model)
                except:
                    if not row.pdb_id == 'XXXX':
                        self.logger.info('The database does not have the model type for %s' % row.pdb_id)

                pdb_chain_to_sym_op[str(row.pdb_id)+"_"+str(row.chain)].add(row.sym_op)


        return pdb_id_to_models, pdb_chain_to_sym_op


    def data(self, pdb, **kwargs):
        """
        Look up all the existing pdbs to process.  Ignores the input pdb variable.
        """

        # map PDB id to models
        # self.logger.info("Starting model query")
        # pdb_id_to_models = self.model_query()

        # map PDB id and chain to symmetry operators
        # self.logger.info("Starting symmetry operator query")
        # pdb_chain_to_sym_op = self.sym_op()

        self.logger.info("Starting model and symmetry operator query")
        pdb_id_to_models, pdb_chain_to_sym_op = self.model_sym_op_query()

        # get pdb ids, chains, molecule type, resolution, technique
        self.logger.info("Starting pdb id, chain, molecule type query")
        with self.session() as session:
            query = session.query(mod.ChainInfo.pdb_id,
                            mod.ChainInfo.chain_name,
                            mod.ChainInfo.entity_macromolecule_type,
                            mod.PdbInfo.resolution,
                            mod.PdbInfo.experimental_technique).\
                            join(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id)
                            # join(mod.UnitInfo, and_(mod.ChainInfo.pdb_id == mod.UnitInfo.pdb_id, mod.ChainInfo.chain_name == mod.UnitInfo.chain))
        result = defaultdict(dict)
        dna_long = ['Polydeoxyribonucleotide (DNA)','polydeoxyribonucleotide']
        rna_long = ['Polyribonucleotide (RNA)','polyribonucleotide']
        hybrid_long = ['polydeoxyribonucleotide/polyribonucleotide hybrid','DNA/RNA Hybrid','NA-hybrid']
        entity_type_check = dna_long + rna_long + hybrid_long + ['Polypeptide(L)','Peptide nucleic acid']

        for row in query:
            if not result.get(row.pdb_id):
                result[row.pdb_id]['chains'] = {}
            result[row.pdb_id]['resolution'] = row.resolution
            result[row.pdb_id]['method'] = row.experimental_technique

            result[row.pdb_id]['model'] = sorted(pdb_id_to_models[row.pdb_id])
            if result[row.pdb_id].get('symmetry'):
                result[row.pdb_id]['symmetry'][row.chain_name] = list(pdb_chain_to_sym_op[str(row.pdb_id)+"_"+str(row.chain_name)])
            else:
                result[row.pdb_id]['symmetry'] = {}
                result[row.pdb_id]['symmetry'][row.chain_name] = list(pdb_chain_to_sym_op[str(row.pdb_id)+"_"+str(row.chain_name)])

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
        """
        Write the data to NA_datafile.pickle

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        copy_file = os.path.join("data",'NA_datafile.pickle') #change name?

        pinfo = self.data(pdb)

        with open(copy_file, 'wb') as fh:
            self.logger.info("process: filename open: %s" % copy_file)
            pickle.dump(pinfo, fh, 2)

        # if the file already exists, do not generate a new file
        # filename = os.path.join("/var/www/html",'NA_datafile.pickle')

        # os.system("ln -s %s %s" % (copy_file, filename))


        pass