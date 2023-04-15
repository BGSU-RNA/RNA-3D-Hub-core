"""
Temporary stage to add modified nucleotides to loop_positions table

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
from sqlalchemy import desc
from sqlalchemy import asc
import re
from pymotifs.units.info import Loader as UnitInfoLoader
from pymotifs.loops.extractor import Loader as InfoLoader
from pymotifs.loops.positions import Loader as PositionLoader
from pymotifs.loops.quality import Loader as QALoader

class Loader(core.Loader):
    merge_data = True
    mark = False
    allow_no_data = True

    dependencies =  set([UnitInfoLoader, InfoLoader, PositionLoader, QALoader]) 


    

    def to_process(self, pdbs, **kwargs):
        '''
        Every pbd is writen twice because data method has two stages and they are processed twice in data method
        
        '''
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).\
                join(mod.LoopPositions,
                     mod.LoopPositions.loop_id == mod.LoopInfo.loop_id).\
                filter(mod.LoopPositions.position_2023 == None).\
                distinct()

            to_use = [(r.pdb_id, i) for r in query for i in [1, 2]]
        
        if len(to_use) > 1:
            return to_use
        else:
            raise core.Skip("No new loops found")


        # test
        # dn_process =  [(pdbs, i) for i in [1, 2]]
        # return dn_process
        # return ['7XD3', '7XD4', '7XD5', '7XD6', '7XD7', '7XNX', '7XNY', '7YG8', '7YG9', '7YGA', '7YGB', '7YGC', '7YGD', '7Z3N', '7Z3O', '7ZTA', '8BQD', '8BQX', '8BR3', '8BSJ', '8BTD', '8BTR', '8BUU', '8C8X', '8C8Y', '8C8Z', '8C90', '8C91', '8C92', '8C93', '8C94', '8C95', '8C96', '8C97', '8C98', '8C99', '8C9A', '8C9B', '8C9C', '8CRX', '8E5T', '8FIZ', '8FOM', '8FON', '8G6W', '8G6X', '8G6Y', '8GHU', '8HD6', '8HD7', '8HFR', '8I7N']

    def has_data(self, pdb_dict, **kwargs):
        """
        This method can query the database after to_process but before data method
        to see if there is already data in the database, and if so, it returns True,
        and the data method does not have to work on this item.
        """

        return False

    def remove(self, pdb, **kwargs):

        """
        with self.session() as session:
            query = session.query(mod.LoopInfo).filter_by(pdb_id=pdb)
            ids = [result.loop_id for result in query]

        if not ids:
            return True

        with self.session() as session:
            return session.query(mod.LoopPositions).\
                filter(mod.LoopPositions.loop_id.in_(ids)).\
                delete(synchronize_session=False)
        """

        return True

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


    def data(self, pdb, **kwargs):
        """
        This method gets called on each item in the list returned by the
        to_process method, one after another.
        We will get one pdb identifier at a time.
        """
        
        stage_1 = pdb[1]
        current_pdb = pdb[0]

        # 1. Loop over all pdbs and adding modified nucleotides to its loop
        if stage_1 == 1:
            print('-------------Working on stage 1-------------')
            with self.session() as session:
                for loop_type in ['HL_%', 'IL_%', 'J3_%']:
                    print('-------------This is', loop_type, '---------------', 'in', current_pdb)
                    if session.query(mod.LoopPositions.loop_id).filter(mod.LoopPositions.loop_id.like(loop_type)).filter(mod.LoopPositions.unit_id.like(current_pdb+'%')).count() == 0:
                        print('No this kind of loops in this pdb')
                    else:
                        subquery = session.query(
                            mod.LoopPositions.loop_positions_id,
                            mod.LoopPositions.loop_id,
                            mod.LoopPositions.position,
                            mod.LoopPositions.border,
                            mod.LoopPositions.unit_id,
                            mod.LoopPositions.position_2023).\
                            filter(mod.LoopPositions.loop_id.like(loop_type)).\
                            filter(mod.LoopPositions.unit_id.like(current_pdb+'%')).\
                            order_by(mod.LoopPositions.loop_id, mod.LoopPositions.position).\
                            subquery()
                        
                        query = session.query(
                            mod.UnitInfo.unit_id,
                            mod.UnitInfo.unit,
                            mod.UnitInfo.model,
                            mod.UnitInfo.sym_op,
                            mod.UnitInfo.chain,
                            mod.UnitInfo.chain_index,
                            mod.UnitInfo.unit_type_id,
                            subquery.c.loop_positions_id,
                            subquery.c.loop_id,
                            subquery.c.position,
                            subquery.c.border,
                            subquery.c.position_2023).\
                            outerjoin(subquery,
                                    subquery.c.unit_id == mod.UnitInfo.unit_id).\
                            filter(mod.UnitInfo.pdb_id == current_pdb).\
                            order_by(mod.UnitInfo.model, mod.UnitInfo.sym_op, mod.UnitInfo.chain, mod.UnitInfo.chain_index)


                        
                        count = 0
                        sanity_check = ["A", "C", "G", "U", "DA", "DC", "DG","DT"]
                        rna_or_dna = ["rna", "dna"]
                        tem_p = 1000 

                        keys_update_existing_nts = ['loop_positions_id', 'position_2023']
                        keys_update_modified_nts = ['loop_positions_id', 'loop_id', 'position', 'bulge', 'flanking', 'border', 'unit_id', 'position_2023']
                        loop_id_start_or_end = {}

                        for r in query:
                            if r.position != None:
                                row = [int(r.loop_positions_id), r.position]
                                entry = dict(zip(keys_update_existing_nts, row))
                                yield mod.LoopPositions(**entry)
                            if r.border == 1:
                                count += 1
                                loop_id = r.loop_id
                                if not loop_id in loop_id_start_or_end:
                                    loop_id_start_or_end[loop_id] = 1
                                else:
                                    loop_id_start_or_end[loop_id] += 1
                                if loop_id_start_or_end[loop_id] % 2 == 1:
                                    loop_id_current = loop_id
                            elif (count % 2 == 1):
                                if r.loop_id == None and r.unit not in sanity_check and r.unit_type_id in rna_or_dna:
                                    tem_p += 1
                                    row = [None, loop_id_current, tem_p, None, 0, 0, r.unit_id, None]
                                    entry = dict(zip(keys_update_modified_nts, row))
                                    yield mod.LoopPositions(**entry)
        # First stage ends

        #2. Order loops that have modified nucleotides
        else:
            print('-------------Working on stage 2-------------')
            with self.session() as session:
                query = session.query(mod.UnitInfo.unit_id,
                                        mod.UnitInfo.chain_index,
                                        mod.UnitInfo.chain,
                                        mod.UnitInfo.unit_type_id,
                                        mod.UnitInfo.unit,
                                        mod.LoopPositions.loop_positions_id,
                                        mod.LoopPositions.loop_id,
                                        mod.LoopPositions.border,
                                        mod.LoopPositions.bulge,
                                        mod.LoopPositions.position_2023,
                                        mod.LoopPositions.flanking,
                                        mod.LoopPositions.position).\
                    join(mod.LoopPositions,mod.LoopPositions.unit_id == mod.UnitInfo.unit_id).\
                    filter(mod.UnitInfo.pdb_id == current_pdb).\
                    order_by(mod.UnitInfo.model, mod.UnitInfo.sym_op, mod.UnitInfo.chain, mod.UnitInfo.chain_index)

            loops_have_modified_nts = set()

            for r in query:
                if r.position > 800:
                    loop_id = r.loop_id
                    loops_have_modified_nts.add(loop_id)
            
            keys_update_existing_nts = ['loop_positions_id', 'position_2023']

            loop_id_to_position_2023 = {}
            if len(loops_have_modified_nts) > 0:
                for loop in loops_have_modified_nts:
                    self.logger.info("Loop %s has found modified nucleotides", loop)
                    print('working on update ' + loop + ' with modified nucleotides, updating their position_2023')
                    with self.session() as session:
                        query = session.query(mod.UnitInfo.unit_id,
                                                mod.UnitInfo.chain_index,
                                                mod.UnitInfo.chain,
                                                mod.UnitInfo.unit_type_id,
                                                mod.UnitInfo.unit,
                                                mod.LoopPositions.loop_positions_id,
                                                mod.LoopPositions.loop_id,
                                                mod.LoopPositions.border,
                                                mod.LoopPositions.bulge,
                                                mod.LoopPositions.position_2023,
                                                mod.LoopPositions.flanking,
                                                mod.LoopPositions.position).\
                            join(mod.LoopPositions,mod.LoopPositions.unit_id == mod.UnitInfo.unit_id).\
                            filter(mod.LoopPositions.loop_id == loop).\
                            order_by(mod.UnitInfo.model, mod.UnitInfo.sym_op, mod.UnitInfo.chain, mod.UnitInfo.chain_index, mod.LoopPositions.loop_id)

                        for r in query:
                            loop_id = r.loop_id
                            if not loop_id in loop_id_to_position_2023:
                                loop_id_to_position_2023[loop_id] = 1
                            else:
                                loop_id_to_position_2023[loop_id] +=1    
                            row = [int(r.loop_positions_id), loop_id_to_position_2023[loop_id]]
                            entry = dict(zip(keys_update_existing_nts, row))
                            yield mod.LoopPositions(**entry)