"""
Program for extracting the information for each atom from a cif file.
"""

import os
import sys
import logging

from models import session
from models import AtomData
from models import PdbAnalysisStatus
from models import PdbUnitIdCorrespondence

from MotifAtlasBaseClass import MotifAtlasBaseClass

from fr3d.cif.reader import CIF


class AtomDataLoader(MotifAtlasBaseClass):
    """A class to load the atom level information from a cif file.
    """

    def __init__(self, commit_size=1000):
        MotifAtlasBaseClass.__init__(self)
        self.commit_size = commit_size
        self.pdb_files_folder = os.path.join(
            self.config['locations']['fr3d_root'], 'FR3D', 'PDBFiles')

    def atoms(self, pdb):
        """Get all Atoms in a pdb. Loads the cif file for the pdb using fr3d
        and then iterates over each atom and yields a AtomData object.

        :pdb: The pdb id to load.
        :yields: AtomData object for each atom.
        """

        cif_file = os.path.join(self.pdb_files_folder, pdb + '.cif')
        reader = CIF(cif_file)
        for structure in reader.structures():
            for residue in structure.residues():
                for atom in residue.atoms():
                    yield AtomData(name=atom.name, x=atom.x, y=atom.y,
                                   z=atom.z, nt_id=atom.component_unit_id())

    def store_pdb(self, pdb):
        """Store all atoms in one pdb.

        :pdb: Pdb ID to store.
        """

        for index, atom in enumerate(self.atoms(pdb)):
            try:
                session.add(atom)
            except:
                logging.warning("Could not add atom: %s")
                return None

            if index % self.commit_size == 0:
                self.__store__(pdb, index)
        self.__store(pdb, index)

    def __store__(self, pdb, index):
        try:
            session.commit()
        except:
            logging.error("Could not store atoms for %s, %s", pdb,
                          index)

    def __remove_atoms__(self, pdbs):
        logging.info("Removing old atom data")

        query = session.query(AtomData).\
            join(PdbUnitIdCorrespondence,
                 PdbUnitIdCorrespondence.unit_id == AtomData.unit_id).\
            filter(PdbUnitIdCorrespondence.pdb.in_(pdbs))
        query.delete(synchronize_session='fetch')

        query = session.query(PdbAnalysisStatus).\
            filter(PdbAnalysisStatus.id.in_(pdbs))
        for status in query:
            logging.info("Resetting status for %s", status.pdb)
            status.atoms = None
            session.merge(status)
        session.commit()

    def __strip_loaded__(self, requested):
        query = session.query(PdbAnalysisStatus).\
            filter(PdbAnalysisStatus.pdb.in_(requested))
        known = set([status.pdb for status in query])

        pdbs = []
        for pdb in requested:
            if pdb in known:
                logging.info("Skipping %s", pdb)
            else:
                pdbs.append(pdb)
        return pdbs

    def __call__(self, pdbs, recalculate=False):
        """Load data for the given atoms. If we have already loaded the given
        pdb then it is skipped.

        :pdbs: List of pdbs to load.
        :recalculate: A flag to force recalculation of the given pdbs. Defaults
        to false.
        """

        recalculate = recalculate or self.config['recalculate']['atoms']
        if recalculate:
            self.__remove_atoms__(pdbs)
        else:
            pdbs = self.__strip_loaded__(pdbs)

        for pdb in pdbs:
            logging.info("Loading PDB: %s", pdb)
            try:
                self.store_pdb(pdb)
            except:
                logging.error("Failed to load: %s", pdb)


def main(pdbs):
    """ Run the loader.

    :pdbs: List of pdbs to use. If None is given then all are imported.
    """
    loader = AtomDataLoader()
    loader.start_logging()

    if not pdbs:
        from PdbInfoLoader import PdbInfoLoader
        P = PdbInfoLoader()
        P.get_all_rna_pdbs()
        pdbs = P.pdbs

    try:
        loader(pdbs, recalculate=True)
    except:
        e = sys.exc_info()[1]
        loader.set_email_subject('Atom loading failed')
        loader._crash(e)

    loader.set_email_subject('Atom loading succesful')
    loader.send_report()

if __name__ == "__main__":
    main(sys.argv[1:])
