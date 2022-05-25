"""Load all unit centers into the database.
"""

import pymotifs.core as core

from pymotifs import models as mod

from pymotifs.units.info import Loader as InfoLoader

from pymotifs.skip_files import SKIP

from sqlalchemy import and_


class Loader(core.SimpleLoader):
    """A class to load all glycosidic centers into the database.
    """

    dependencies = set([InfoLoader])
    allow_no_data = True
    
    def to_process(self, pdbs, **kwargs):
        a = '5A9Z 1E7X 4V4G 204D 1ZFR 226D 229D 199D 1QSX 1JU0 1JUU 1M6A 1JRW 1MSF 1MSE 5Z9W 175D 1DSI 1BCE 1BCB 6GV4 1BJH 1BJD 1KSE 1XUE 230D 1AP1 1EVO 1UQE 1UQD 1UQG 1UQF 1UQA 1UQC 1UQB 2STT 2STW 1BBX 1HZ0 1HZ2 2ARG 1FZL 1FZX 1FZS 3GAT 1I7V 1BPS 1K9L 1K9H 1EVM 1EVN 1CK8 1AC9 1LA8 1AC7 1LAS 1LAE 1LAI 2NEO 1QDF 1QDK 1QDH 1QDI 6QNQ 1B3P 1FYY 1FYI 1FKY 1FKZ 1BWG 1BWT 177D 1OKA 1JO1 1EW1 139D 1IG4 1EEG 1EEK 1G1N 1KBD 2KBD 1FJB 5MJV 4Z92 1FJ5 170D 149D 1QE7 4GAT 2GAT 1KVH 1I5V 1IV6 2DAU 1IDX 1RCS 1KXS 1GCC 1G22 1ZHU 132D 1LCD 1LCC 1SJK 1HM1 1CQO 1AGK 1WAN 1BUT 1BUF 1BUB 6NUT 148D 1HRY 1HRZ 1IEK 1IEY 1KKW 1KKV 1EU6 1CN8 1CN9 2NVS 214D 5GAT 1A66 1ELN 1FV8 1CR3 1QCH 1D6D 1T42 1SKP 1I34 1B0S 103D 1EZN 1D42 1AO9 1AO1 1D69 4WZD 2NYO 1HO6 1SLS 1HOQ 1CS2 1DK9 1DK6 1C0Y 1N14 1D70 1AL9 1YUI 1YUJ 1G7Z 1DHH 1QSK 141D 1DXA 1BDZ 1C11 1AMD 169D 105D 1DJD 5M74 1MXK 1D68 1GIZ 1GIP 202D 1G80 104D 1DUF 1JS7 1MYQ 2EZE 2EZD 2EZG 2EZF 1C32 6I2N 1C95 1EOR 1GJ2 1KB1 1AXL 1NLE 1AXV 1AXU 1AXP 186D 1DRN 1AX6 1AX7 1GTC 1DZS 1JRV 203D 1B4Y 1BJ6 1C38 1C35 1C34 1OLD 5MSF 1HT7 1HT4 1SAA 1JS5 107D 1BX5 1K1H 1H8J 1K1R 1EL2 1F5E 3REC 5A79 5A7A 1EMQ 5APO 1GKW 1GKV 1FC8 1G14 1NK2 1NK3 1NYZ 1HE0 1HE6 1DSM 1HWV 106D 1GN7 179D 1DGO 1FQP 4KBD 1DAU 1K2J 1K2K 140D 1K29 1TNE 1J5N 108D 1J5K 1DL4 1AFF 1D3X 1AFZ 1ECU 1CFL 1LVS 1F4S 1CX3 1AT4 1BHR 2HDC 201D 1BE5 1L0R 1J46 6MSF 1HVO 1HVN 1D20 1EU2 1AF1 171D 1F3S 5FN1 1D18 1D19 1KBM 1LWA 3KBD 1AUL 1AG3 1CYZ 4WR6 6H5Q 1AGU 1LEJ 1AGZ 4WRO 1G4D 1AU6 1AGH 1A8N 1QBY 1A8W 1A6B 1A6H 1GJ0 1K4X 1A83 1N1K 135D 1N17 142D 1II1 1E6T 1ESS 1COC 1CK5 1AXO 225D 1G5L 1G5E 1G5D 1AGO 7MSF 1C7U 1HDW 1JVC 1JVE 1BN9 1TAN 1HX4 5AA0 143D 1JDG 136D 1LAQ 5AFI 1DB6 1K5F 1JJN 1JJM 156D 1JJP 1KR8 1K5E'
        
        return a.split()




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

        return list(pdb_ids - existed_ids - set(SKIP))
        # return ['4V3P']

    def query(self, session, pdb):
        return session.query(mod.UnitCenters).\
            filter(mod.UnitCenters.pdb_id == pdb).\
            filter(mod.UnitCenters.name == 'glycosidic')

    def data(self, pdb, **kwargs):
        structure = self.structure(pdb)
        for residue in structure.residues():
            for name in residue.centers.definitions():
                center = residue.centers[name]
                if len(center) == 3 and name == 'glycosidic':
                    yield mod.UnitCenters(unit_id=residue.unit_id(),
                                          name=name,
                                          pdb_id=pdb,
                                          x=float(center[0]),
                                          y=float(center[1]),
                                          z=float(center[2]))
