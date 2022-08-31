"""Stage to populate the units.info table.

This module contains a loader to load all unit level information into the
database.
"""

import itertools as it

import pymotifs.core as core

from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs import utils
# from pymotifs.download import Downloader
# from pymotifs.pdbs.info import Loader as PdbLoader
# from pymotifs.units import Loader as UnitInfoLoder
from sqlalchemy import and_
from collections import defaultdict
from Bio.Alphabet import ThreeLetterProtein

AA = [seq.upper() for seq in ThreeLetterProtein().letters]

class Loader(core.Loader):
    merge_data = True
    mark = False


    dependencies = set([])
    """The dependencies for this stage."""

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,mod.UnitInfo.pdb_id).\
                filter(mod.UnitInfo.unit_type_id == None)
            if not query.count():
                raise core.Skip("Skipping summary, no new correspondences")
        d = defaultdict(dict)
        for result in query:
            if d.get(result.pdb_id):
                d[result.pdb_id].append(result.unit_id)
            else:
                d[result.pdb_id] = []
                d[result.pdb_id].append(result.unit_id)
        return list(d.values())

    def remove(self, pdb_id, **kwargs):
        self.logger.info("Not removing anything, recompute all correspondence")

    def has_data_2(self, unit_id, **kwargs):

        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter(mod.UnitInfo.unit_id == unit_id).\
                filter(mod.UnitInfo.unit_type_id != None)
            return bool(query.count())
    def has_data(self, pdb_dict, **kwargs):
        return 0

        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter(mod.UnitInfo.unit_id.in_(pdb_dict.values())).\
                filter(mod.UnitInfo.unit_type_id != None)
            return bool(query.count())

    def current(self, unit_id): # current row info
        """Get the current data for the correspondence.
        """

        with self.session() as session:
            info = session.query(mod.UnitInfo).get(unit_id)
            return utils.row2dict(info)


    def type_query(self, unit_ids, **kwargs):
        d = {}
        structure = self.structure(unit_ids[0][:4])
        # print(units.component_type(self.structure(unit_ids[:4]).unit_id()))
        # print(units.component_type(self.structure(unit_ids[:4]).residues('5Z9W|1|A|ALA|297||||P_35')))
        for base in structure.residues(polymeric=None):
            d.update({base.unit_id():units.component_type(base)})
        for unit_id in unit_ids:
            if unit_id.split('|')[3] in AA:
                d.update({unit_id:'aa'})
        return d


    def data(self, units_list, **kwargs):
        # print(units_list)
        type_query_dict = self.type_query(units_list)
        # print(type_query_dict)
        for unit_id in units_list:
            data = self.current(unit_id)
            try:
                row_update = type_query_dict[unit_id]
                data['unit_type_id'] = row_update
            # print(data)
            except:
                self.logger.info("The pipeline does not have unit_type_id for %s"%unit_id)
            yield mod.UnitInfo(**data)

