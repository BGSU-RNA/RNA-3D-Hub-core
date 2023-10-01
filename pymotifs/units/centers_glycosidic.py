"""Load all unit centers into the database.
"""

import pymotifs.core as core

from pymotifs import models as mod
import itertools as it
from pymotifs.units.info import Loader as InfoLoader

from pymotifs.skip_files import SKIP
from pymotifs.utils import units
from sqlalchemy import and_


class Loader(core.SimpleLoader):
    """
    A class to load all glycosidic centers into the database.
    Ignores pdbs passed in.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True

    def to_process(self, pdbs, **kwargs):
        # return ["1UVK"]

        with self.session() as session:
            query_all_pdb_ids = session.query(mod.ChainInfo.pdb_id).distinct()
            query_glycosidic_present = session.query(mod.UnitCenters.pdb_id).\
                                        filter(mod.UnitCenters.name == 'glycosidic').distinct()
                                        #filter(and_(mod.UnitCenters.name == 'base',mod.UnitCenters.name != 'glycosidic')).distinct()

        pdb_ids = set()
        for result in query_all_pdb_ids:
            pdb_ids.add(result.pdb_id)
        existing_ids = set()
        for result in query_glycosidic_present:
            existing_ids.add(result.pdb_id)

        if len(pdb_ids - existing_ids - set(SKIP)) == 0:
            raise core.Skip("No new glycosidic centers")

        self.logger.info('%s'%list(pdb_ids - existing_ids - set(SKIP)))
        return sorted(pdb_ids - existing_ids - set(SKIP))
        # return ['4V3P']

    def query(self, session, pdb):
        return session.query(mod.UnitCenters).\
            filter(mod.UnitCenters.pdb_id == pdb).\
            filter(mod.UnitCenters.name == 'glycosidic')

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


    def data(self, pdb, **kwargs):

        # with self.session() as session:
        #     existed_unit_ids = session.query(mod.UnitInfo.unit_id).\
        #                         filter(mod.UnitInfo.pdb_id == pdb)
        # unit_ids_list = []
        # for result in existed_unit_ids:
        #     unit_ids_list.append(result.unit_id)
        # print(unit_ids_list)

        structure = self.structure(pdb)
        for residue in structure.residues():
            # if residue.unit_id() not in unit_ids_list:
            #     print(residue.unit_id())
            #     print(residue.pdb)
            #     print(residue.model)
            #     print(residue.chain)
            #     print(residue.sequence)
            #     print(residue.number)
            #     print(getattr(residue, 'alt_id', None))
            #     print(residue.insertion_code)
            #     print(residue.symmetry)
            #     print(residue.index)
            #     print(self.type(residue))
            #     yield mod.UnitInfo(unit_id=residue.unit_id(),
            #             pdb_id=residue.pdb,
            #             model=residue.model,
            #             chain=residue.chain,
            #             unit=residue.sequence,
            #             number=residue.number,
            #             alt_id=getattr(residue, 'alt_id', None),
            #             ins_code=residue.insertion_code,
            #             sym_op=residue.symmetry,
            #             chain_index=residue.index,
            #             unit_type_id=self.type(residue))
            for name in residue.centers.definitions():
                center = residue.centers[name]
                if len(center) == 3 and name == 'glycosidic':
                    yield mod.UnitCenters(unit_id=residue.unit_id(),
                                          name=name,
                                          pdb_id=pdb,
                                          x=float(center[0]),
                                          y=float(center[1]),
                                          z=float(center[2]))
