from __future__ import with_statement
"""
Compute the ordering of units in a PDB file.

Example:
    python unit_ordering_loader.py
"""

import os
import logging

from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import session, PdbUnitOrdering, PdbModifiedCorrespondecies

from Bio.PDB.PDBParser import PDBParser

logger = logging.getLogger(__name__)


class UnitOrderingLoader(MotifAtlasBaseClass):

    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.commit_every = 100
        self.file_types = ['.pdb', '.pdb1']
        root = self.config['locations']['fr3d_root']
        self.pdb_files_folder = os.path.join(root, 'FR3D', 'PDBFiles')

        self.known = ['A', 'C', 'G', 'U']
        query = session.query(PdbModifiedCorrespondecies.modified_unit)
        for entry in query.all():
            self.known.append(entry.modified_unit)

    def import_ordering(self, pdbs, recalculate=False):
        """Imports the unit id ordering for the given pdb files. Will only
        import things if their unit is in the list of known units.
        """
        logger.info('Importing ordering')
        analyze = pdbs
        if recalculate or self.config['recalculate']['ordering']:
            self.__delete_ordering__(pdbs)
        else:
            analyze = self.filter_out_analyzed_pdbs(analyze, 'unit_ordering')

        for pdb in analyze:
            for file_type in self.file_types:
                try:
                    self.__import_file__(pdb, file_type)
                except:
                    logger.error("Failed ordering import: %s%s",
                                 pdb, file_type)

    def __import_file__(self, pdb, extension):
        """Attempt to import a file. If the file does not exists, a warning is
        logged otherwise importing occurs.
        """
        filename = os.path.join(self.pdb_files_folder, pdb + extension)

        if not os.path.exists(filename):
            logger.info("Skipping missing file: %s", filename)
            return None

        try:
            with open(filename, 'r') as raw:
                pdb_type = 'AU'
                if extension[-1] != 'b':
                    pdb_type = 'BA' + extension[-1]
                ordering = self.__pdb_ordering__(raw, pdb, pdb_type)
                self.__store__(ordering)
                self.mark_pdb_as_analyzed(pdb, 'unit_ordering')
        except:
            logger.critical('Crash on: %s', filename)

    def __store__(self, data):
        """Store the results of generating an ordering.
        """
        logger.info("Storing ordered ids")
        try:
            for count, (unit_id, info) in enumerate(data.items()):
                try:
                    entry = PdbUnitOrdering(nt_id=unit_id, pdb=info['pdb'],
                                            index=info['index'])
                    session.add(entry)
                    if count % self.commit_every == 0:
                        session.commit()
                except:
                    logger.error("Failed to add: %s", unit_id)
            session.commit()
            logger.info("ID ordering added")
        except:
            logger.info("Could not commit all ids")

    def __pdb_ordering__(self, raw, pdb_id, pdb_type):
        """Generate a dict of the form: { unit_id: {index: index, pdb: pdb }
        for all nucleotides in the given structure. Nucleotides are identified
        by being in the list of known units in self.known.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, raw)
        data = {}
        index = 0
        for model in structure:
            model_id = model.get_id() + 1
            for chain in model:
                chain_id = chain.get_id()
                for residue in chain:
                    name = residue.resname.strip()
                    if name in self.known:
                        res_id = residue.get_id()
                        id_data = [structure.get_id(), pdb_type, model_id,
                                   chain_id, res_id[1], name, res_id[2]]
                        id_data = [str(part).strip() for part in id_data]
                        unit_id = '_'.join(id_data)
                        data[unit_id] = {'index': index, 'pdb': pdb_id}
                        index += 1

        return data

    def __delete_ordering__(self, pdbs):
        logger.info("Deleting ordering")

        entries = session.query(PdbUnitOrdering). \
            filter(PdbUnitOrdering.pdb.in_(pdbs))

        entries.delete(synchronize_session='fetch')

if __name__ == "__main__":
    import sys
    pdbs = sys.argv[1:]
    loader = UnitOrderingLoader()
    loader.start_logging()
    loader.load(pdbs)
