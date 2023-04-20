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
    """A class to load all glycosidic centers into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True
    
    def to_process(self, pdbs, **kwargs):
        # new dna structures
        # return ['1EHL']
        # return ['4WT1', '4WU1', '4WZO', '4X62', '4X64', '4X65', '4X66', '4XEJ', '4Y4O', '4Y4P', '4YHH', '4YPB', '4YY3', '4YZV', '4Z3Q', '4Z3R', '4Z3S', '4ZER', '4ZSN', '5A9Z', '5BR8', '5CZP', '5D8B', '5DFE', '5DOX', '5DOY', '5E7K', '5E81', '5EL4', '5EL5', '5EL6', '5EL7', '5F8K', '5FDU', '5FDV', '5FN1', '5IB7', '5IB8', '5IBB', '5IMQ', '5IMR', '5IWA', '5J30', '5J3C', '5J4B', '5J4C', '5J4D', '5LMN', '5LMO', '5LMP', '5LMQ', '5LMR', '5LMS', '5LMT', '5LMU', '5LMV', '5LZB', '5LZC', '5LZD', '5LZE', '5LZF', '5OMW', '5ON2', '5ON3', '5ONH', '5OT7', '5UQ7', '5UQ8', '5V8I', '5VP2', '5VPO', '5VPP', '5W4K', '5WIS', '5WIT', '5WNP', '5WNQ', '5WNR', '5WNS', '5WNT', '5WNU', '5WNV', '5ZLU', '6B4V', '6BOH', '6BOK', '6BUW', '6BZ6', '6BZ7', '6BZ8', '6C5L', '6CAE', '6CAO', '6CAP', '6CAQ', '6CAR', '6CAS', '6CFJ', '6CFK', '6CFL', '6FKR', '6GSJ', '6GSK', '6GSL', '6HIV', '6HIX', '6MKN', '6MPI', '6N1D', '6N9E', '6N9F', '6ND5', '6ND6', '6NDK', '6NY6', '6O97', '6OF1', '6OF6', '6OJ2', '6OPE', '6ORD', '6OTR', '6OXA', '6OXI', '6Q95', '6QNQ', '6QNR', '6UCQ', '6UO1']
        # return ['1BJ6', '1CGM', '1FJF', '1FJG', '1GIY', '1H3E', '1HAO', '1HAP', '1HNW', '1HNX', '1HNZ', '1HR0', '1HUT', '1IBK', '1IBL', '1IBM', '1J5E', '1ML5', '1N32', '1N33', '1N34', '1N36', '1PFI', '1RMV', '1S1H', '1VS9', '1VSA', '1VSP', '1VTM', '1VVJ', '1VVM', '1VVO', '1VVQ', '1VVS', '1VVU', '1VVW', '1VVY', '1VW0', '1VX9', '1VXJ', '1VXL', '1VXN', '1VXQ', '1VXT', '1VY1', '1VY3', '1VY4', '1VY5', '1VY6', '1VY7', '1XMO', '1XMQ', '1XNQ', '1XNR', '1YL3', '1YL4', '2B64', '2B66', '2B9M', '2B9N', '2B9O', '2B9P', '2E5L', '2F4V', '2HGJ', '2HGQ', '2HGU', '2J00', '2J01', '2J02', '2J03', '2JL5', '2JL6', '2JL7', '2JL8', '2NR0', '2NZ4', '2OM3', '2TMV', '2UU9', '2UUA', '2UUB', '2UUC', '2UXB', '2UXC', '2V46', '2V47', '2V48', '2V49', '2VQE', '2VQF', '2WDG', '2WDH', '2WDI', '2WDJ', '2WDK', '2WDL', '2WDM', '2WDN', '2WH1', '2WH2', '2WH3', '2WH4', '2WRI', '2WRJ', '2WRK', '2WRL', '2WRN', '2WRO', '2WRQ', '2WRR', '2WYY', '2X9R', '2X9S', '2X9T', '2X9U', '2XEA', '2XFZ', '2XG0', '2XG1', '2XG2', '2XQD', '2XQE', '2XSY', '2XTG', '2XUX', '2XUY', '2Y0U', '2Y0V', '2Y0W', '2Y0X', '2Y0Y', '2Y0Z', '2Y10', '2Y11', '2Y12', '2Y13', '2Y14', '2Y15', '2Y16', '2Y17', '2Y18', '2Y19', '2ZM6', '2ZZM', '3A3A', '3ADB', '3ADC', '3ADD', '3D5B', '3D5D', '3DD2', '3F1F', '3F1H', '3FIC', '3FIN', '3G8S', '3G8T', '3G95', '3G96', '3G9C', '3HL2', '3HUW', '3HUX', '3HUY', '3HUZ', '3I8F', '3I8I', '3I9C', '3I9E', '3J06', '3KIQ', '3KIR', '3KIS', '3KIT', '3KIU', '3KIW', '3KIX', '3KIY', '3KNH', '3KNI', '3KNJ', '3KNK', '3KNL', '3KNM', '3KNN', '3KNO', '3L3C', '3MRZ', '3MS1', '3OGE', '3OGY', '3OH5', '3OH7', '3OHC', '3OHD', '3OHJ', '3OHK', '3OHY', '3OHZ', '3OI0', '3OI1', '3OI2', '3OI3', '3OI4', '3OI5', '3OTO', '3PDM', '3PYO', '3PYR', '3PYT', '3PYV', '3RG5', '3TVE', '3TVH', '3UXQ', '3UXR', '3UXS', '3UXT', '3UYE', '3UYG', '3UZ1', '3UZ2', '3UZ8', '3UZ9', '3UZF', '3UZK', '3UZN', '3V22', '3V23', '3V24', '3V25', '3V26', '3V27', '3V28', '3V29', '3V2C', '3V2D', '3V2E', '3V2F', '3V6W', '3V6X', '3W1K', '3W3S', '3ZGZ', '3ZJT', '3ZJU', '3ZJV', '3ZN7', '3ZN9', '3ZND', '3ZNE', '3ZVO', '3ZVP', '4ABR', '4ABS', '4AQ7', '4AQY', '4ARC', '4ARI', '4AS1', '4B3M', '4B3R', '4B3S', '4B3T', '4B8F', '4B8G', '4B8H', '4B8I', '4BTC', '4BTD', '4BYB', '4BYC', '4BYD', '4BYE', '4CQN', '4CR1', '4DH9', '4DHA', '4DHB', '4DHC', '4DR1', '4DR2', '4DR3', '4DR4', '4DR5', '4DR6', '4DR7', '4DUY', '4DUZ', '4DV0', '4DV1', '4DV2', '4DV3', '4DV4', '4DV5', '4DV6', '4DV7', '4EJ9', '4EJA', '4EJB', '4EJC', '4G5L', '4G5N', '4G5U', '4G5W', '4JI0', '4JI1', '4JI2', '4JI3', '4JI4', '4JI5', '4JI6', '4JI7', '4JI8', '4JUW', '4JUX', '4JV5', '4K0L', '4KBT', '4KBU', '4KBV', '4KBW', '4KCY', '4KCZ', '4KD0', '4KD2', '4KD8', '4KD9', '4KDA', '4KDB', '4KDG', '4KDH', '4KDJ', '4KDK', '4KFH', '4KFI', '4KFK', '4KFL', '4KVB', '4KX0', '4KX2', '4L47', '4L6J', '4L6L', '4L71', '4LEL', '4LF4', '4LF5', '4LF6', '4LF7', '4LF8', '4LF9', '4LFA', '4LFB', '4LFC', '4LFZ', '4LNT', '4LSK', '4LT8', '4NIA', '4NVU', '4NVW', '4NVX', '4NVY', '4NW0', '4NW1', '4NXM', '4NXN', '4OX9', '4P6F', '4P70', '4QCM', '4QCO', '4QCP', '4QCQ', '4QCS', '4QCT', '4QCU', '4QCW', '4QCX', '4QCY', '4QD0', '4QD1', '4QJS', '4QJT', '4QS0', '4QS1', '4QS2', '4QS3', '4RB5', '4RB7', '4RB8', '4RB9', '4RBB', '4RBC', '4RBD', '4RBF', '4RBG', '4RBH', '4RBJ', '4RBK', '4RQE', '4RQF', '4TUA', '4TUB', '4TUC', '4TUD', '4TUE', '4UDV', '4V42', '4V4B', '4V4I', '4V4J', '4V4P', '4V4R', '4V4S', '4V4T', '4V4X', '4V4Y', '4V4Z', '4V51', '4V5A', '4V5C', '4V5D', '4V5E', '4V5F', '4V5G', '4V5J', '4V5K', '4V5L', '4V5M', '4V5N', '4V5P', '4V5Q', '4V5R', '4V5S', '4V63', '4V67', '4V68', '4V6A', '4V6F', '4V6G', '4V7J', '4V7K', '4V7L', '4V7M', '4V7P', '4V7W', '4V7X', '4V7Y', '4V7Z', '4V83', '4V84', '4V87', '4V8A', '4V8B', '4V8C', '4V8D', '4V8E', '4V8F', '4V8G', '4V8H', '4V8I', '4V8J', '4V8N', '4V8O', '4V8Q', '4V8U', '4V8X', '4V90', '4V95', '4V97', '4V9A', '4V9B', '4V9H', '4V9I', '4V9J', '4V9K', '4V9L', '4V9M', '4V9N', '4V9Q', '4V9R', '4V9S', '4W29', '4W2B', '4W2D', '4W2E', '4W2F', '4W2G', '4W2H', '4W2I', '4W4G', '4WPO', '4WQ1', '4WQF', '4WQR', '4WQU', '4WQY', '4WRA', '4WSD', '4WSM', '4WT1', '4WU1', '4WZO', '4X62', '4X64', '4X65', '4X66', '4XEJ', '4Y4O', '4Y4P', '4YHH', '4YPB', '4YY3', '4YZV', '4Z3Q', '4Z3R', '4Z3S', '4ZER', '4ZSN', '5A9Z', '5BR8', '5CZP', '5D8B', '5DFE', '5DOX', '5DOY', '5E7K', '5E81', '5EL4', '5EL5', '5EL6', '5EL7', '5F8K', '5FDU', '5FDV', '5FN1', '5IB7', '5IB8', '5IBB', '5IMQ', '5IMR', '5IWA', '5J30', '5J3C', '5J4B', '5J4C', '5J4D', '5LMN', '5LMO', '5LMP', '5LMQ', '5LMR', '5LMS', '5LMT', '5LMU', '5LMV', '5LZB', '5LZC', '5LZD', '5LZE', '5LZF', '5OMW', '5ON2', '5ON3', '5ONH', '5OT7', '5UQ7', '5UQ8', '5V8I', '5VP2', '5VPO', '5VPP', '5W4K', '5WIS', '5WIT', '5WNP', '5WNQ', '5WNR', '5WNS', '5WNT', '5WNU', '5WNV', '5ZLU', '6B4V', '6BOH', '6BOK', '6BUW', '6BZ6', '6BZ7', '6BZ8', '6C5L', '6CAE', '6CAO', '6CAP', '6CAQ', '6CAR', '6CAS', '6CFJ', '6CFK', '6CFL', '6FKR', '6GSJ', '6GSK', '6GSL', '6HIV', '6HIX', '6MKN', '6MPI', '6N1D', '6N9E', '6N9F', '6ND5', '6ND6', '6NDK', '6NY6', '6O97', '6OF1', '6OF6', '6OJ2', '6OPE', '6ORD', '6OTR', '6OXA', '6OXI', '6Q95', '6QNQ', '6QNR', '6UCQ', '6UO1']
        # return ['1VS9', '1FJF']
        # return ['5GAT','4KBD','4GAT','3REC','3KBD','3GAT','2STW','2STT','2NYO','2NVS']
        # return ['1EO4', '6QNQ']
        # return ['5GAT']
        # return ['7MSF', '1EO4', '6MSF', '6NUT', '6QNQ', '5MJV', '5MSF', '4Z92', '5FN1', '5M74']
        # return ['7MSF', '1EO4', '6MSF', '6NUT', '1VS9', '2I1C', '6QNQ', '1FJF', '5MJV', '5MSF', '4Z92', '5FN1', '5M74']
        # return ['6I2N', '4WR6', '6H5Q', '4WRO', '4WZD', '7MSF', '1EO4', '6MSF', '6NUT', '5A79', '5A7A', '5APO', '1VS9', '2I1C', '6QNQ', '5AFI', '1FJF', '5AA0', '5MJV', '5MSF', '5Z9W', '4Z92', '5FN1', '6GV4', '5M74']
        # return ['1EO4']
        # return ['135D','136D','1BJ6','1D6D','1IDX','1II1','1JDG','1NLE','1QE7','1ZFR','3HL2','4ZDO','4ZDP']



        with self.session() as session:
            query_all_pdb_ids = session.query(mod.ChainInfo.pdb_id).distinct()
            query_glycosidic_present = session.query(mod.UnitCenters.pdb_id).\
                                        filter(mod.UnitCenters.name == 'glycosidic').distinct()
                                        #filter(and_(mod.UnitCenters.name == 'base',mod.UnitCenters.name != 'glycosidic')).distinct()

        # pdb_ids = []
        # for id in query_all_pdb_ids:
        #     pdb_ids.append(str(id)[2:6])
        # existed_ids = []
        # for id in query_glycosidic_present:
        #     existed_ids.append(str(id)[2:6])
        # passing_pdb_ids = set()
        # for pdb_id in pdb_ids:
        #     if (pdb_id not in SKIP) and (pdb_id not in existed_ids):
        #         passing_pdb_ids.add(pdb_id)

        pdb_ids = set()
        for id in query_all_pdb_ids:
            pdb_ids.add(str(id)[2:6])
        existed_ids = set()
        for id in query_glycosidic_present:
            existed_ids.add(str(id)[2:6])        

        if len(pdb_ids - existed_ids - set(SKIP)) == 0:
            raise core.Skip("No new glycosidic centers")


        # return [pdb_id for pdb_id in result for skip_id in SKIP if pdb_id != skip_id]
        # print(len([pdb_id for pdb_id in result for skip_id in SKIP if pdb_id != skip_id]))
        # print(len(pdb_ids))
        # print(len(passing_pdb_ids))
        # print(len(SKIP))
        # self.logger.info('length of pdb_ids:%s length of processing:%s length of existed ids:%s' % (len(pdb_ids),len(passing_pdb_ids),len(existed_ids)))
        self.logger.info('%s'%list(pdb_ids - existed_ids - set(SKIP)))
        return list(pdb_ids - existed_ids - set(SKIP))
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
            ##### do not delete the following lines because they can fill in rows for unit_info table.
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
