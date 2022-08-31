"""Stage to populate the units.info table.

This module contains a loader to load all unit level information into the
database.
"""

import itertools as it

import pymotifs.core as core

from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs.download import Downloader
from pymotifs.pdbs.info import Loader as PdbLoader
from sqlalchemy import and_


class Loader(core.SimpleLoader):
    """The loader that will populate the unit_info table in the database.
    """

    dependencies = set([Downloader, PdbLoader])
    """The dependencies for this stage."""

    def query(self, session, pdb):
        """Create a query for all units for the given PDB.

        Parameters
        ----------
        session : Session
            The session to use.
        pdb : str
            The PDB ID to use.

        Returns
        -------
        query : sqlalchemy.orm.query.Query
            A query to find all units from the given structure.
        """
        return session.query(mod.UnitInfo).filter_by(pdb_id=pdb)

    def type(self, unit):
        """Compute the component type, ie A, C, G, U is RNA, DA, DC, etc is DNA
        and so forth.

        Parameters
        ----------
        unit : Component
            The unit to get the component for

        Returns
        -------
        component_type : str
            The component type.
        """
        return units.component_type(unit)

    # def to_process(self, pdbs, **kwargs):
    #     try:
    #         for i in pdbs:
    #             structure = self.structure(i).residues(polymeric=None)
    #             j = 1
    #             for base in structure:
    #                 print(base.unit_id(),base.type,base.sequence,base.symmetry, self.type(base),base.chain, self.unit_type_None_fixed(base), 'row',j)
    #                 j = j+1
    #     except:
    #         print('running')
    #     return(pdbs)
    def unit_type_None_fixed(self, nt):
        """we may not need to check the entity_macromolecule_type in the chain info table
        we can get the RNA linking or DNA linking information by nt.type
        we added more if condictions to the units.component_type function, so this function is no longer needed.
        """
        if nt.type == 'RNA linking':
            return 'rna'
        elif nt.type == 'DNA linking':
            return 'dna'
        else:
            return None


    def as_unit(self, nt):
        """Turn a `Component` into a `UnitInfo`.

        Parameters
        ----------
        nt : Component
            The `Component` to turn into a `UnitInfo`.

        Returns
        -------
        unit : UnitInfo
            The `Component` as a `UnitInfo`
        """
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
        """Compute the data to store. This will extract all components from the
        structure, include water, ligands and other non-polymers and create
        `UnitInfo` objects.

        Parameters
        ----------
        pdb : str
            The PDB id to compute the units for.

        Returns
        -------
        data : iterator
            An iterable over all units in the structure.
        """

        structure = self.structure(pdb)

        return it.imap(self.as_unit, structure.residues(polymeric=None))
