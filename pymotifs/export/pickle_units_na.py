"""
Module for export of unit center/rotation data in pickle format for FR3D
"""

import numpy as np
import os
import pickle

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import DATA_FILE_DIRECTORY

from pymotifs.chains.info import Loader as ChainLoader
from pymotifs.units.center_rotation import Loader as CenterRotationsLoader
from pymotifs.exp_seq.mapping import Loader as MappingLoader
from pymotifs.exp_seq.positions import Loader as PositionLoader
from pymotifs.ife.info import Loader as IfeInfoLoader


class Exporter(core.Loader):
    """
    Export unit data in pickle format, one file per chain.
    """

    # General Setup
    compressed = False
    mark = False

    # Write out all files, not just ones for current pdb files
    write_all_files = True
    write_all_files = False

    # consider all files, not just current RNA PDB files
    consider_all_files = True
    consider_all_files = False

    dependencies = set([ChainLoader, CenterRotationsLoader,
                        PositionLoader, IfeInfoLoader, MappingLoader])

    def remove():
        pass


    def filename(self, chain_tuple, **kwargs):
        """
        Create the filename for the given chain.

        Parameters
        ----------
        chain_tuple : tuple
            The components (PDB ID, model, and chain) of the
            chain name (hyphen-deliminted concatenation)
            for which to create a file.

        Returns
        -------
        filename : str
            The path to write to.
        """

        # TO DO: put the important directories into the config

        pdb = chain_tuple[0]
        mdl = chain_tuple[1]
        chn = chain_tuple[2]

        chain_string = pdb + '-' + str(mdl) + '-' + chn

        return os.path.join(DATA_FILE_DIRECTORY,'units',chain_string + "_NA.pickle")


    def to_process(self, pdbs, **kwargs):
        """
        Look up the list of chains to process.
        Only process pdbs in the given list.
        If you want to produce *all* such files, use self.write_all_files.

        Parameters
        ----------
        pdbs : list

        Returns
        -------
        (pdb_id, model, chain) : tuple
            The components of the chains to be processed.
        """

        with self.session() as session:
            query = session.query(
                       mod.UnitInfo.pdb_id,
                       mod.UnitInfo.model,
                       mod.UnitInfo.chain).\
                   filter(mod.UnitInfo.unit_type_id.in_(['rna','dna'])).\
                   distinct()

            if not self.write_all_files and not self.consider_all_files:
                query = query.filter(mod.UnitInfo.pdb_id.in_(pdbs))

            na_chains = set([(r.pdb_id, str(r.model), r.chain) for r in query])
            self.logger.info('Found %d nucleic acid chains in unit_info' % len(na_chains))

        na_type_list = ['Polydeoxyribonucleotide (DNA)','Polyribonucleotide (RNA)','DNA/RNA Hybrid','polyribonucleotide','polydeoxyribonucleotide','polydeoxyribonucleotide/polyribonucleotide hybrid']

        with self.session() as session:
            query = session.query(
                       mod.ChainInfo.pdb_id,
                       mod.ChainInfo.chain_name).\
                   filter(mod.ChainInfo.entity_macromolecule_type.in_(na_type_list)).\
                   distinct()

            if not self.write_all_files and not self.consider_all_files:
                query = query.filter(mod.ChainInfo.pdb_id.in_(pdbs))

            na_chains_table = set([(r.pdb_id, "1", r.chain_name) for r in query])
            self.logger.info('Found %d nucleic acid chains in chain_info' % len(na_chains_table))

        na_chains |= na_chains_table
        self.logger.info('Found %d nucleic acid chains in all' % len(na_chains))

        # get a list of all files in directory that end with _NA.pickle
        directory = os.path.join(DATA_FILE_DIRECTORY,'units')
        files = [f for f in os.listdir(directory) if f.endswith('_NA.pickle')]
        written_chains = set()
        for f in files:
            written_chains.add(tuple(f.split("_")[0].split('-')))

        self.logger.info('Found %d na chains already written' % len(written_chains))

        if self.write_all_files:
            need_to_write = na_chains
        else:
            need_to_write = sorted(na_chains - written_chains)

        if len(need_to_write) == 0:
            raise core.Skip("No new pdb files to process")

        return sorted(need_to_write)


    def has_data(self, entry, *args, **kwargs):

        if self.write_all_files:
            # literally write all files
            return False

        #self.logger.info("has_data: entry: %s" % str(entry))
        filename = self.filename(entry)
        #self.logger.info("has_data: filename: %s" % filename)

        # temporary to overwrite some/all files
        if os.path.exists(filename):
            self.logger.info("has_data: filename %s exists" % filename)
            return True
        else:
            self.logger.info("has_data: filename %s is missing" % filename)
            return False


    def data(self, chain_tuple, **kwargs):
        """
        Get all unit listings for the given chain, centers and
        rotations, and format them for convenient use by FR3D.

        Parameters
        ----------
        chain_tuple : tuple
            The chain for which to look up centers and rotation data.

        Returns
        -------
        resultset : list of dicts
        """

        pdb = chain_tuple[0]
        mdl = chain_tuple[1]
        chn = chain_tuple[2]

        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,
                               mod.UnitInfo.chain_index,
                               mod.UnitCenters.x,
                               mod.UnitCenters.y,
                               mod.UnitCenters.z,
                               mod.UnitRotations.cell_0_0,
                               mod.UnitRotations.cell_0_1,
                               mod.UnitRotations.cell_0_2,
                               mod.UnitRotations.cell_1_0,
                               mod.UnitRotations.cell_1_1,
                               mod.UnitRotations.cell_1_2,
                               mod.UnitRotations.cell_2_0,
                               mod.UnitRotations.cell_2_1,
                               mod.UnitRotations.cell_2_2).\
                     distinct().\
                     join(mod.UnitCenters, mod.UnitInfo.unit_id == mod.UnitCenters.unit_id).\
                     join(mod.UnitRotations, mod.UnitInfo.unit_id == mod.UnitRotations.unit_id).\
                     filter(mod.UnitCenters.name == 'glycosidic').\
                     filter(mod.UnitInfo.pdb_id == pdb).\
                     filter(mod.UnitInfo.model == mdl).\
                     filter(mod.UnitInfo.chain == chn).\
                     order_by(mod.UnitInfo.chain_index)

            # removed this filter on 2025-02-28 because glycosidic center is enough
            #  filter(mod.UnitInfo.unit_type_id.in_(['rna','dna'])).\

            # removed this filter on 2025-02-27 for solitary nucleotides like 2JEF|1|A|DGT|1343
            #  filter(mod.UnitInfo.chain_index != None).\

            self.logger.debug("cenrot: query built")

            units = []
            order = []
            cntrs = []
            rttns = []
            rsset = []

            for row in query:
                units.append(row.unit_id)
                if row.chain_index == None:
                    # solitary nucleotides that are not covalently linked to a chain
                    order.append(None)
                else:
                    order.append(row.chain_index - 1)
                cntrs.append(np.asarray([row.x, row.y, row.z]))
                rttns.append(np.asarray([np.asarray([row.cell_0_0, row.cell_0_1, row.cell_0_2]),
                             np.asarray([row.cell_1_0, row.cell_1_1, row.cell_1_2]),
                             np.asarray([row.cell_2_0, row.cell_2_1, row.cell_2_2])]))

            rsset = [ units, order, cntrs, rttns ]

            self.logger.debug("cenrot: units: %s" % str(units))
            self.logger.debug("cenrot: order: %s" % str(order))
            self.logger.debug("cenrot: cntrs: %s" % str(cntrs))
            self.logger.debug("cenrot: rttns: %s" % str(rttns))
            self.logger.debug("cenrot: rsset: %s" % str(rsset))

            return rsset


    def process(self, entry, **kwargs):
        """
        Load centers/rotations data for the given chain.

        Parameters
        ----------
        entry : object
            The entry to process.
        **kwargs : dict
            Generic keyword arguments.
        """
        # if len(self.data(entry)[0]) != 0:
        # write the file even if the list is empty, to assert that no nucleotides have centers
        # 9Z80 is an example, where all nucleotides are UNK.  9Q3G has all N.
        filename = self.filename(entry)

        uinfo = self.data(entry)

        with open(filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(uinfo, fh, 2)

