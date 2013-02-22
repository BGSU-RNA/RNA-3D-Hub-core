#from future import with_statment

"""
Compute the ordering of units in a PDB file.

Example:
    python unit_ordering_loader.py
"""

import os
import logging

from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import session, PdbUnitOrdering

from Bio.PDB.PDBParser import PDBParser


class UnitOrderingLoader(MotifAtlasBaseClass):
    commit_every = 100
    known = []
    file_types = ['.pdb', '.pdb1']

    def load(self, pdbs, recalculate=False):
        """Imports the unit id ordering for the given pdb files. Will only
        import things if their unit is in the list of known units.
        """
        logging.info('Importing ordering')
        analyze = pdbs
        if recalculate or self.config['recalculate']['ordering']:
            self.__delete_ordering__(pdbs)
        else:
            analyze = self.filter_out_analyzed_pdbs(analyze, 'unit_ordering')

        for pdb in analyze:
            path = os.path.join(self.pdb_files_folder, pdb)
            for file_type in self.file_types:
                try:
                    self.__import_file__(path, file_type)
                except:
                    logging.error("Failed ordering import: %s%s",
                                  pdb, file_type)

    def __import_file__(self, pdb, extension):
        """Attempt to import a file. If the file does not exists, a warning is
        logged otherwise importing occurs.
        """
        filename = pdb + extension

        if not os.path.exists(filename):
            logging.info("Skipping missing file: %s", filename)

        try:
            with open(filename, 'r') as raw:
                ordering = self.__pdb_ordering__(raw, pdb, extension)
                self.__store__(ordering)
        except:
            logging.cricical('Crash on: %s', filename)

    def __store__(self, data):
        """Store the results of generating an ordering.
        """
        logging.info("Storing ordered ids")
        for count, (unit_id, index) in enumerate(data.entries()):
            entry = PdbUnitOrdering(unit_id=unit_id, index=index)
            try:
                session.add(entry)
                if count % self.commit_every == 0:
                    session.commit()
            except:
                logging.error("Failed to add: %s", unit_id)
                pass
        session.commit()
        logging.info("ID ordering added")

    def __pdb_ordering__(self, raw, pdb_id, pdb_type):
        """Generate a dict of the form: { unit_id: index } for all nucleotides
        in the given structure. Nucleotides are identified by being in the list
        of known units.
        """
        parser = PDBParser()
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
                        data[unit_id] = index
                        index += 1

        return data

    def __delete_ordering__(self, pdbs):
        logging.info("Deleting ordering")

        entries = session.query(PdbUnitOrdering). \
            filter(PdbUnitOrdering.pdb.in_(pdbs))

        entries.delete(synchronize_session='fetch')

if __name__ == "__main__":
    import sys
    pdbs = sys.argv[1:]
    loader = UnitOrderingLoader()
    loader.start_logging()
    loader.load(pdbs)
