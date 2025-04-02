"""
Find structures that do not have glycosidic centers calculated, and load them.
"""

import itertools as it

import pymotifs.core as core
from pymotifs import models as mod
from pymotifs.units.info import Loader as InfoLoader
from pymotifs.skip_files import SKIP
from pymotifs.utils import units

class Loader(core.SimpleLoader):
    """
    A class to load all glycosidic centers into the database.
    Ignores pdbs passed in.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True
    mark = False

    def to_process(self, pdbs, **kwargs):
        """
        If pdbs is a list of just one pdb id, process that one.
        Otherwise, find all pdbs that do not have glycosidic
        """

        if len(pdbs) == 1:
            return pdbs
        else:
            with self.session() as session:
                query_all_pdb_ids = session.query(mod.ChainInfo.pdb_id).distinct()
                pdb_ids = set()
                for result in query_all_pdb_ids:
                    pdb_ids.add(result.pdb_id)
            self.logger.info("Total number of pdb files: %s" % len(pdb_ids))

            # this query takes a long time
            with self.session() as session:
                query_glycosidic_present = session.query(mod.UnitCenters.pdb_id).\
                                            filter(mod.UnitCenters.name == 'glycosidic').distinct()

                pdb_with_glycosidic = set()
                for result in query_glycosidic_present:
                    pdb_with_glycosidic.add(result.pdb_id)
            self.logger.info("Number of pdb files with glycosidic centers: %s" % len(pdb_with_glycosidic))

            if len(pdb_ids - pdb_with_glycosidic - set(SKIP)) == 0:
                raise core.Skip("All pdb files have glycosidic centers")

            self.logger.info('%s' % list(pdb_ids - pdb_with_glycosidic - set(SKIP)))

            return sorted(pdb_ids - pdb_with_glycosidic - set(SKIP))


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

        # load the structure
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
                    self.logger.info('Adding glycosidic center for %s' % residue.unit_id())
                    yield mod.UnitCenters(unit_id=residue.unit_id(),
                                          name=name,
                                          pdb_id=pdb,
                                          x=float(center[0]),
                                          y=float(center[1]),
                                          z=float(center[2]))
