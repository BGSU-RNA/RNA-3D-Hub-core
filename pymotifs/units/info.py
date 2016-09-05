import itertools as it

import pymotifs.core as core

from Bio.Alphabet import ThreeLetterProtein

from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs.download import Downloader
from pymotifs.pdbs.info import Loader as PdbLoader

AA = [seq.upper() for seq in ThreeLetterProtein().letters]


class Loader(core.SimpleLoader):
    update_gap = False
    dependencies = set([Downloader, PdbLoader])

    def query(self, session, pdb):
        """Create a query for all units for the given PDB.

        :param Session session: The session to use.
        :param str pdb: The PDB Id to use.
        :returns: A query to find all units from the given structure.
        """
        return session.query(mod.UnitInfo).filter_by(pdb_id=pdb)

    def type(self, unit):
        """Compute the component type, ie A, C, G, U is RNA, DA, DC, etc is DNA
        and so forth.

        :param Component unit: The unit to get the component for
        :returns: The component type.
        """
        return units.component_type(unit)

    def as_unit(self, nt):
        return mod.UnitInfo(unit_id=nt.unit_id(),
                            pdb_id=nt.pdb,
                            model=nt.model,
                            chain=nt.chain,
                            unit=nt.sequence,
                            number=nt.number,
                            alt_id=getattr(nt, 'alt_id', None),
                            ins_code=nt.insertion_code,
                            sym_op=nt.symmetry,
                            chain_index=nt.index,
                            unit_type_id=self.type(nt))

    def data(self, pdb, **kwargs):
        structure = self.structure(pdb)
        return it.imap(self.as_unit, structure.residues(polymeric=None))
