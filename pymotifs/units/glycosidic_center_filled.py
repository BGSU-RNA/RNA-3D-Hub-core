"""Load all unit centers into the database.
"""

from collections import defaultdict
import pymotifs.core as core

from pymotifs import models as mod
import itertools as it
from pymotifs.units.info import Loader as InfoLoader

from pymotifs.skip_files import SKIP
from pymotifs.utils import units
from sqlalchemy import and_


class Loader(core.SimpleLoader):
    """A class to load all glycosidic centers into the database.
    """

    dependencies = set([])
    allow_no_data = True
    
    def to_process(self, pdbs, **kwargs):
        # return ['1Q96']
        print(111)
        with self.session() as session:
            query_glycosidic_fill = session.query(mod.UnitCenters.pdb_id, mod.UnitCenters.unit_id, mod.UnitCenters.name).\
                filter(mod.UnitCenters.pdb_id == '1Q96')

        unit_id_with_base = set()
        unit_id_with_glycosidic = set()
        pdb_id_to_unit_ids  = defaultdict(dict)
        for row in query_glycosidic_fill:
            if row.name == 'base':
                unit_id_with_base.add(row.unit_id)
            if row.name == 'glycosidc':
                unit_id_with_glycosidic.add(row.unit_id)
        unit_id_without_glycosidic = unit_id_with_base - unit_id_with_glycosidic

        for unit_id in unit_id_without_glycosidic:
            if pdb_id_to_unit_ids.get(str(unit_id).split('|')[0]):
                pdb_id_to_unit_ids[str(unit_id).split('|')[0]].append(unit_id)
            else:
                pdb_id_to_unit_ids[str(unit_id).split('|')[0]] = []
                pdb_id_to_unit_ids[str(unit_id).split('|')[0]].append(unit_id)



        # extend_pdb_id = set()
        # for unit_id in unit_id_without_glycosidic:
        #     extend_pdb_id.add(str(unit_id).split("|")[0])
        print(222)
        print(pdb_id_to_unit_ids.keys())
        return pdb_id_to_unit_ids



    def query(self, session, pdb):
        return session.query(mod.UnitCenters).\
            filter(mod.UnitCenters.pdb_id == pdb).\
            filter(mod.UnitCenters.name == 'glycosidic')

    # def type(self, unit):
    #     """Compute the component type, ie A, C, G, U is RNA, DA, DC, etc is DNA
    #     and so forth.

    #     Parameters
    #     ----------
    #     unit : Component
    #         The unit to get the component for

    #     Returns
    #     -------
    #     component_type : str
    #         The component type.
    #     """
    #     print(444)
    #     return units.component_type(unit)      


    def data(self, as_dict, **kwargs):
        print(333)
        print(as_dict.keys()[0])

        pdb = list(as_dict.keys())[0]
        structure = self.structure(pdb)


        structure = self.structure(pdb)
        for residue in structure.residues():

            for name in residue.centers.definitions():
                center = residue.centers[name]
                if len(center) == 3 and name == 'glycosidic':
                    if residue.unit_id() not in as_dict.values():
                        yield mod.UnitCenters(unit_id=residue.unit_id(),
                                            name=name,
                                            pdb_id=pdb,
                                            x=float(center[0]),
                                            y=float(center[1]),
                                            z=float(center[2]))
