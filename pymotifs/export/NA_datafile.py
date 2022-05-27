"""Module for export of resolution and methods of existing pdb files
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
        """Create the filename for the given PDB.

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

        self.logger.info("filename: filename: %s" % filename)
        ## change /html to /unit
        return os.path.join("/var/www/units",filename)
        # return os.path.join("/usr/local/pipeline/hub-core/logs",filename)

    def to_process(self, pdbs, **kwargs):
        """Ignore the pdbs input.
            The return value is just a list value.
            We do this because we only want to run this script one time even if we have a lot of pdbs.
        """
        return ['1']

    def model_query(self, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitInfo.pdb_id,mod.UnitInfo.chain,mod.UnitInfo.model).distinct()
                            # filter(and_(mod.UnitInfo.pdb_id == pdb, mod.UnitInfo.chain == chain_name)).first()
                            # filter(mod.UnitInfo.pdb_id == pdb).distinct()
        # return str(query)[1]
        # the output of result will be like this format, (12L)
        # thus, (, ), and L need to be removed.
        # return str(query).replace('(','').replace(')','').replace('L','').replace(',','')
        result = defaultdict(dict)
        for row in query:
            # if result[row.pdb_id].get(row.chain):
            #     result[row.pdb_id][row.chain].append(row.model)
            # else:
            #     result[row.pdb_id][row.chain] = []
            #     result[row.pdb_id][row.chain].append(row.model)
            result[row.pdb_id][row.chain]=str(row.model).replace('L','')
        return result      
        

    def data(self, pdb, **kwargs):
        """
            Look up all the existed pdbs to process.  Ignores the pdb input.
        """
        model_type = self.model_query()

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
        hybrid_long = ['polydeoxyribonucleotide/polyribonucleotide hybrid','DNA/RNA Hybrid']
        entity_type_check = dna_long + rna_long + hybrid_long + ['Polypeptide(L)','Peptide nucleic acid']

        for row in query:
            if not result.get(row.pdb_id):
                result[row.pdb_id]['chains'] = {}
            result[row.pdb_id]['resolution'] = row.resolution
            result[row.pdb_id]['method'] = row.experimental_technique
            # try:
            #     if result[row.pdb_id].get('model'):
            #         result[row.pdb_id]['model'].append(model_type[row.pdb_id][row.chain_name])
            #     else:
            #         result[row.pdb_id]['model'] = []
            #         result[row.pdb_id]['model'].append(model_type[row.pdb_id][row.chain_name])
            # except:
            #     self.logger.info('The database currently do not have the model type for %s with chain %s' % (row.pdb_id,row.chain_name))


            if row.entity_macromolecule_type == 'Polypeptide(L)':
                if result[row.pdb_id]['chains'].get('protein'):
                    result[row.pdb_id]['chains']['protein'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['protein'] = []
                    result[row.pdb_id]['chains']['protein'].append(row.chain_name)
            ###########
                try:
                    if result[row.pdb_id].get('protein_model'):
                        result[row.pdb_id]['protein_model'].append(model_type[row.pdb_id][row.chain_name])
                    else:
                        result[row.pdb_id]['protein_model'] = []
                        result[row.pdb_id]['protein_model'].append(model_type[row.pdb_id][row.chain_name])
                except:
                    self.logger.info('The database currently do not have the protein model type for %s with chain %s' % (row.pdb_id,row.chain_name))      
            elif (row.entity_macromolecule_type in dna_long):
                if result[row.pdb_id]['chains'].get('DNA'):
                    result[row.pdb_id]['chains']['DNA'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['DNA'] = []
                    result[row.pdb_id]['chains']['DNA'].append(row.chain_name)
            ###########
                try:
                    if result[row.pdb_id].get('DNA_model'):
                        result[row.pdb_id]['DNA_model'].append(model_type[row.pdb_id][row.chain_name])
                    else:
                        result[row.pdb_id]['DNA_model'] = []
                        result[row.pdb_id]['DNA_model'].append(model_type[row.pdb_id][row.chain_name])
                except:
                    self.logger.info('The database currently do not have the DNA model type for %s with chain %s' % (row.pdb_id,row.chain_name))                       
            elif (row.entity_macromolecule_type in rna_long):
                if result[row.pdb_id]['chains'].get('RNA'):
                    result[row.pdb_id]['chains']['RNA'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['RNA'] = []
                    result[row.pdb_id]['chains']['RNA'].append(row.chain_name)
            ###########
                try:
                    if result[row.pdb_id].get('RNA_model'):
                        result[row.pdb_id]['RNA_model'].append(model_type[row.pdb_id][row.chain_name])
                    else:
                        result[row.pdb_id]['RNA_model'] = []
                        result[row.pdb_id]['RNA_model'].append(model_type[row.pdb_id][row.chain_name])
                except:
                    self.logger.info('The database currently do not have the RNA model type for %s with chain %s' % (row.pdb_id,row.chain_name))   
            elif (row.entity_macromolecule_type in hybrid_long):
                if result[row.pdb_id]['chains'].get('hybrid'):
                    result[row.pdb_id]['chains']['hybrid'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['hybrid'] = []
                    result[row.pdb_id]['chains']['hybrid'].append(row.chain_name)
            ###########
                try:
                    if result[row.pdb_id].get('hybrid_model'):
                        result[row.pdb_id]['hybrid_model'].append(model_type[row.pdb_id][row.chain_name])
                    else:
                        result[row.pdb_id]['hybrid_model'] = []
                        result[row.pdb_id]['hybrid_model'].append(model_type[row.pdb_id][row.chain_name])
                except:
                    self.logger.info('The database currently do not have the hybrid model type for %s with chain %s' % (row.pdb_id,row.chain_name))   
            elif (row.entity_macromolecule_type == 'Peptide nucleic acid'):
                if result[row.pdb_id]['chains'].get('PNA'):
                    result[row.pdb_id]['chains']['PNA'].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains']['PNA'] = []
                    result[row.pdb_id]['chains']['PNA'].append(row.chain_name)
            ###########
                try:
                    if result[row.pdb_id].get('PNA_model'):
                        result[row.pdb_id]['PNA_model'].append(model_type[row.pdb_id][row.chain_name])
                    else:
                        result[row.pdb_id]['PNA_model'] = []
                        result[row.pdb_id]['PNA_model'].append(model_type[row.pdb_id][row.chain_name])
                except:
                    self.logger.info('The database currently do not have the PNA model type for %s with chain %s' % (row.pdb_id,row.chain_name))  
            elif (row.entity_macromolecule_type not in entity_type_check):
                if result[row.pdb_id]['chains'].get(row.entity_macromolecule_type):
                    result[row.pdb_id]['chains'][row.entity_macromolecule_type].append(row.chain_name)
                else:
                    result[row.pdb_id]['chains'][row.entity_macromolecule_type] = []
                    result[row.pdb_id]['chains'][row.entity_macromolecule_type].append(row.chain_name)
            ###########
                other_type_model = str(row.entity_macromolecule_type) + '_model'
                try:
                    if result[row.pdb_id].get(other_type_model):
                        result[row.pdb_id][other_type_model].append(model_type[row.pdb_id][row.chain_name])
                    else:
                        result[row.pdb_id][other_type_model] = []
                        result[row.pdb_id][other_type_model].append(model_type[row.pdb_id][row.chain_name])
                except:
                    self.logger.info('The database currently do not have the %s model type for %s with chain %s' % (row.entity_macromolecule_type,row.pdb_id,row.chain_name)) 
            


            

        print(result)

        return result



    def process(self, pdb, **kwargs):
        """Load centers/rotations data for the given IFE-chain.

        Parameters
        ----------
        **kwargs : dict
            Generic keyword arguments.
        """

        webroot = self.config['locations']['fr3d_pickle_base'] + "/"

        filename = self.filename()
        copy_file = os.path.join("pickle-FR3D",'NA_datafile.pickle')

        pinfo = self.data(pdb)


        with open(filename, 'wb') as fh:
            self.logger.info("process: filename open: %s" % filename)
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(pinfo, fh, 2)

        os.system("rsync -u %s %s" % (filename, webroot))

        with open(copy_file, 'wb') as fh:
            self.logger.info("process: filename open: %s" % copy_file)
            pickle.dump(pinfo, fh, 2)

        os.system("rsync -u %s %s" % (copy_file, webroot))


        pass


