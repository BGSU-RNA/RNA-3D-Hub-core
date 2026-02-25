"""
Stage to populate the units.info table.

This module contains a loader to load all unit level information into the
database.
"""

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs.download import Downloader
from pymotifs.pdbs.info import Loader as PdbLoader

class Loader(core.SimpleLoader):
    """
    The loader that will populate the unit_info table in the database.
    """

    # the dependencies for this stage
    dependencies = set([Downloader, PdbLoader])

    allow_no_data = True
    use_marks = False

    # check for units that were missed and fill them in
    fill_in_missing = True
    fill_in_missing = False

    def query(self, session, pdb):
        """
        Create a query for all units for the given PDB.

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

        if self.fill_in_missing:
            # return an empty query so we process all pdb ids
            return session.query(mod.UnitInfo).filter_by(pdb_id='nonexistent_pdb_id')
        else:
            return session.query(mod.UnitInfo).filter_by(pdb_id=pdb)


    def type(self, unit):
        """
        Compute the component type, ie A, C, G, U is RNA, DA, DC, etc is DNA
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


    def as_unit(self, nt, unit_id="", symmetry=""):
        """
        Turn a `Component` into a `UnitInfo`.

        Parameters
        ----------
        nt : Component
            The `Component` to turn into a `UnitInfo`.

        Returns
        -------
        unit : UnitInfo
            The `Component` as a `UnitInfo`
        """
        # self.logger.info("unit_id: %s, pdb: %s, model: %s, chain: %s, unit: %s, number: %s, alt_id: %s, ins_code: %s, sym_op: %s, chain_index: %s, unit_type_id: %s" % (nt.unit_id(), nt.pdb, nt.model, nt.chain, nt.sequence, nt.number, getattr(nt, 'alt_id', None), nt.insertion_code, nt.symmetry, nt.index, self.type(nt)))

        if not unit_id:
            unit_id = nt.unit_id()

        if not symmetry:
            symmetry = nt.symmetry

        return mod.UnitInfo(unit_id=unit_id,
                            pdb_id=nt.pdb,
                            model=nt.model,
                            chain=nt.chain,
                            unit=nt.sequence,
                            number=nt.number,
                            alt_id=getattr(nt, 'alt_id', None),
                            ins_code=nt.insertion_code,
                            sym_op=symmetry,
                            chain_index=nt.index,
                            unit_type_id=self.type(nt))


    def convert_symmetry(self,u,old_prefix,new_prefix):
        # convert ASM_ to P_ for saving in the database
        if old_prefix:
            f = u.split("|")
            if len(f) == 9:
                if f[8].startswith(old_prefix):
                    f[8] = f[8].replace(old_prefix,new_prefix)
                    u = "|".join(f)
        return u


    def data(self, pdb, **kwargs):
        """
        Compute the data to store. This will extract all components from the
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

        old_prefix = ""
        new_prefix = ""

        if self.fill_in_missing:
            # identify unit ids and symmetry operators already present
            with self.session() as session:
                query = session.query(mod.UnitInfo).filter_by(pdb_id=pdb)
                existing = set([u.unit_id for u in query])

                database_symmetry = set([u.sym_op.split("_")[0] for u in query])
                self.logger.info('database_symmetry is %s' % database_symmetry)

                if 'P' in database_symmetry and not 'ASM' in database_symmetry:
                    old_prefix = "ASM_"
                    new_prefix = "P_"
                    self.logger.info('Changing unit ids being added from %s to %s' % (old_prefix,new_prefix))

        else:
            existing = set()

        structure = self.structure(pdb)

        entries = []
        for unit in structure.residues():
            u = unit.unit_id()
            s = unit.symmetry
            if old_prefix:
                u = self.convert_symmetry(u,old_prefix,new_prefix)
                s = s.replace(old_prefix,new_prefix)
            if not u in existing:
                if self.fill_in_missing:
                    self.logger.info("Filling in missing unit %s with symmetry %s" % (u,s))
                entries.append(self.as_unit(unit,u,s))

        self.logger.info("Found %4d new units for %s" % (len(entries), pdb))

        return entries
