"""
Fill in unit_type_id in unit_info table, where possible
"""

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.utils import units
from pymotifs import utils


class Loader(core.Loader):
    merge_data = True
    mark = False

    dependencies = set([])

    def to_process(self, pdbs, **kwargs):

        if len(pdbs) == 1:
            # get all units in pdbs[0] with unit_type_id equal to NULL
            with self.session() as session:
                query = session.query(mod.UnitInfo.unit_id,mod.UnitInfo.unit).\
                    filter(mod.UnitInfo.unit_type_id == None).\
                    filter(mod.UnitInfo.pdb_id == pdbs[0])
        else:
            # get all units with unit_type_id equal to NULL
            with self.session() as session:
                query = session.query(mod.UnitInfo.unit_id,mod.UnitInfo.unit).\
                    filter(mod.UnitInfo.unit_type_id == None)

                if not query.count():
                    raise core.Skip("Skipping summary, no new correspondences")

        # map distinct values of unit to list of unit ids
        unit_to_unit_ids = {}
        for result in query:
            unit = result.unit
            if unit in unit_to_unit_ids:
                unit_to_unit_ids[unit].append(result.unit_id)
            else:
                unit_to_unit_ids[unit] = []
                unit_to_unit_ids[unit].append(result.unit_id)

        all_tuples = list(unit_to_unit_ids.items())

        # sort tuples from longest list to shortest
        all_tuples.sort(key=lambda x: len(x[1]), reverse=True)

        for unit, unit_ids in all_tuples:
            self.logger.info("Processing %s with %d instances" % (unit, len(unit_ids)))

        # return a list of tuples of unit and list of unit ids
        return all_tuples


    def remove(self, pdb_id, **kwargs):
        self.logger.info("Not removing anything, recompute all correspondence")


    def has_data(self, pdb_dict, **kwargs):
        return 0

        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter(mod.UnitInfo.unit_id.in_(pdb_dict.values())).\
                filter(mod.UnitInfo.unit_type_id != None)
            return bool(query.count())


    def current(self, unit_id): # current row info
        """
        Get the current data for the unit id from unit_info table
        """

        with self.session() as session:
            info = session.query(mod.UnitInfo).get(unit_id)
            return utils.row2dict(info)


    def get_type_from_cif(self, unit_id, **kwargs):
        """
        Read a .cif file to determine the unit_type_id for each residue
        Doesn't really help with the difficult cases
        """

        structure = self.structure(unit_id.split("|")[0])

        for base in structure.residues():
            print(base.unit_id(),base.type)
            if base.unit_id() == unit_id:
                return base.type

        return None


    def data(self, unit_and_unit_ids, **kwargs):
        """
        Loop over units and the unit ids where they occur
        Look up the unit type
        When found, update the unit_info table
        """

        unit = unit_and_unit_ids[0]
        unit_ids = unit_and_unit_ids[1]

        if not unit:
            raise core.Skip("No unit found")

        # unit_type_id = self.get_type_from_cif(unit_ids[0])

        # call the function where we focus our effort to get these right
        # use just the name of the unit
        unit_type_id = units.component_type(None, unit)

        if unit_type_id:
            self.logger.info("Setting unit_type_id to %s for %s" % (unit_type_id, unit))
            self.logger.info("unit_ids: %s" % unit_ids)
            for unit_id in unit_ids:
                data = self.current(unit_id)
                data['unit_type_id'] = unit_type_id
                yield mod.UnitInfo(**data)
                # raise core.Skip("Skipping writing to the database during testing")
        else:
            raise core.Skip("No unit_type_id found for unit %s" % unit)
