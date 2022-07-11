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


class Loader(core.Loader):
    merge_data = True
    mark = False


    dependencies = set([])
    """The dependencies for this stage."""

    def to_process(self, pdbs, **kwargs):
        """We transform the list of pdbs into the list of correspondences that
        have not yet been summarized.

        :param list pdb: The list of pdb ids. Currently ignored.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of correspondence ids to process.
        """

        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id).\
                filter(mod.UnitInfo.unit_type_id == None)
            if not query.count():
                raise core.Skip("Skipping summary, no new correspondences")

        return [result.unit_id for result in query]
    def remove(self, unit_id, **kwargs):
        self.logger.info("Not removing anything, recompute all correspondence")

    def has_data(self, unit_id, **kwargs):

        with self.session() as session:
            query = session.query(mod.UnitInfo).\
                filter(mod.UnitInfo.unit_id == unit_id).\
                filter(mod.UnitInfo.unit_type_id != None)
            return bool(query.count())

    def current(self, unit_id):
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
    def type_query(self, unit_id, **kwargs):
        structure = self.structure(unit_id[:4])
        d = {}
        for base in structure.residues():
            d.update({base.unit_id():self.type(base)})
        #print(d)
        return d                 


    def data(self, unit_id, **kwargs):
        """Compute the summary for the given correspondence id. This will
        update the entry with the counts of match, mismatch and such.
        """

        data = self.current(unit_id)
        # print(data)
        # print("111",data['unit_type_id'],data['unit_id'])
        try:
            row_update = self.type_query(unit_id)[unit_id]
            data['unit_type_id'] = row_update
        except:
            self.logger.info("The pipeline do not have unit_type_id for metal ion")

        # if row_update:
        #     data['unit_type_id'] = row_update
        # else:
        #     data['unit_type_id'] = None
        # print("222",data['unit_type_id'],data['unit_id'])
        # print(data)
        # print(row_update)


        return mod.UnitInfo(**data)

