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


class Loader(core.Loader):
    merge_data = True
    mark = False


    dependencies = set([])
    """The dependencies for this stage."""

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id).\
                filter(mod.UnitInfo.unit_type_id == None)
            if not query.count():
                raise core.Skip("Skipping summary, no new correspondences")

        return [result.unit_id for result in query]

    def to_process_new(self, pdbs, **kwargs):
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,mod.UnitInfo.pdb_id).\
                filter(mod.UnitInfo.unit_type_id == None).\
                filter(mod.UnitInfo.pdb_id == '101D')
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

    # def to_process_Z(self, pdbs, **kwargs):

    #     with self.session() as session:
    #         query = session.query(mod.UnitInfo.pdb_id).\
    #             filter(mod.UnitInfo.unit_type_id == None).\
    #             filter(mod.UnitInfo.pdb_id == '101D')
    #         if not query.count():
    #             raise core.Skip("Skipping summary, no new correspondences")

    #     return [result.pdb_id for result in query]


    def remove(self, pdb_id, **kwargs):
        self.logger.info("Not removing anything, recompute all correspondence")

    def has_data(self, unit_id, **kwargs):

        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter(mod.UnitInfo.unit_id == unit_id).\
                filter(mod.UnitInfo.unit_type_id != None)
            return bool(query.count())
    def has_data_new(self, pdb_dict, **kwargs):
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


    def type_query_old(self, unit_id, **kwargs):
        structure = self.structure(unit_id[:4])
        d = {}
        for base in structure.residues():
            d.update({base.unit_id():self.type(base)})
        return d

    def type_query(self, unit_id, **kwargs):
        print(unit_id)
        if unit_id.split('|')[3].upper() in ['CU', 'FE', 'MG', 'NI', 'MN', 'K', 'NA', 'MO', 'CO', 'ZN', 'W', 'CA', 'V']:
            return 'ion'
        base = self.structure(unit_id[:4]).residues(unit_id)
        # if base.unit_id().split('|')[3].upper() in ['CU', 'FE', 'MG', 'NI', 'MN', 'K', 'NA', 'MO', 'CO', 'ZN', 'W', 'CA', 'V']:
        #     return {base.unit_id():'ion'}
        # print(base.unit_id().split('|')[3])
        d = self.type(base)
        return d    

    # def type_query(self, pdb_id, **kwargs):
    #     structure = self.structure(pdb_id)
    #     d = {}
    #     for base in structure.residues():
    #         d.update({base.unit_id():self.type(base)})
    #     return d             

    # def data(self, pdb_ids, **kwargs):
    #     for pdb in pdb_ids:
    #         with self.session() as session:
    #             query = session.query(mod.UnitInfo.unit_id).\
    #                 filter(mod.UnitInfo.unit_type_id == None).\
    #                 filter(mod.UnitInfo.pdb_id == pdb)
    #         type_dict = self.type_query(pdb_id)
    #         for row in query:
    #             data = self.current(row.unit_id)
    #             row_update = type_dict[row.unit_id]
    #             data['unit_type_id'] = row_update


    #         yield mod.UnitInfo(**data)


    def data(self, unit_id, **kwargs):
        """Compute the summary for the given correspondence id. This will
        update the entry with the counts of match, mismatch and such.
        """
        # print(unit_id)
        # for unit_id in unit_ids:
        # data = self.current(unit_id)
        # row_update = self.type_query(unit_id)
        # data['unit_type_id'] = row_update
        data = self.current(unit_id)
        try:
            row_update = self.type_query(unit_id)
            data['unit_type_id'] = row_update
        except:
            self.logger.info("The pipeline do not have unit_type_id for %s"%unit_id)


        return mod.UnitInfo(**data)

    def data_new(self, pdb_dict, **kwargs):
        print(pdb_dict)
        for unit_id in pdb_dict.values():
            data = self.current(unit_id)
            try:
                row_update = self.type_query(unit_id)
                data['unit_type_id'] = row_update
            except:
                self.logger.info("The pipeline do not have unit_type_id for %s"%unit_id)

##############
            yield mod.UnitInfo(**data)

