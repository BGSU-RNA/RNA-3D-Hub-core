"""
Read HL, IL, J3, etc. from a file created by NA_pairwise_annotations.py
Read loops currently stored in the database, if any, and compare
If old loops need to be deprecated, do that by indicating it in loop_info
"""

import glob
import os
from sqlalchemy import or_

from fr3d.modified.mapping import modified_base_to_parent

from pymotifs import core
from pymotifs import models as mod

from pymotifs.pdbs.info import Loader as PdbLoader
from pymotifs.units.info import Loader as UnitLoader
from pymotifs.interactions.loader import Loader as InteractionLoader


class Loader(core.SimpleLoader):
    merge_data = True
    allow_no_data = True
    dependencies = set([PdbLoader, UnitLoader, InteractionLoader])
    save_loops = True
    use_marks = False
    mark = False

    # location of loops from NA_pairwise_annotations
    loop_location = '/usr/local/pipeline/hub-core/data/loops'
    done_location = '/usr/local/pipeline/hub-core/data/loops_done'

    def print_dictionary(self, d):
        for key in sorted(d.keys()):
            if not "_sa_" in key:
                print("%25s %s" % (key, d[key]))

    def dictionary_as_text(self, d):
        t = "\n"
        for key in sorted(d.keys()):
            if not "_sa_" in key and not "flanking" in key and not key == "position":
                t += "%-21s: %s\n" % (key, d[key])
        return t

    def to_process(self, pdbs, **kwargs):

        # development code, to identify pdb ids where it is safe to delete loops
        if False:
            # query the database to find all loop ids in loop_info where the length of
            # the field sort_name is greater than 0
            with self.session() as session:
                query = session.query(mod.LoopInfo.loop_id).\
                    filter(mod.LoopInfo.sort_name != '')
                saved_loop_ids = set([r.loop_id for r in query])

            # query the motif atlas for loop ids that are in the motif atlas
            with self.session() as session:
                query = session.query(mod.MlLoops.loop_id).distinct()
                motif_atlas_loop_ids = set([r.loop_id for r in query])

            print('Found %5d loops in loop_info' % len(saved_loop_ids))
            print('Found %5d loops in motif_atlas' % len(motif_atlas_loop_ids))
            print('Found these %5d loops in both:' % (len(saved_loop_ids & motif_atlas_loop_ids)))
            print(sorted(saved_loop_ids & motif_atlas_loop_ids))

        # for temporary overrides, speeds up testing
        pdbs_needed = ['6TNA']
        pdbs_needed = ['4TNA']
        pdbs_needed = ['4V9F']
        pdbs_needed = ['1FFK']
        pdbs_needed = ['1S72']
        pdbs_needed = ['1VWW']
        pdbs_needed = ['1DK1']
        pdbs_needed = ['1DUH']
        pdbs_needed = ['1F7Y']
        pdbs_needed = []         # leave blank to use all pdbs

        if len(pdbs_needed) == 0:
            # get the list of pdb ids that already have loops in loop_info table
            with self.session() as session:
                query = session.query(mod.LoopInfo.pdb_id).distinct()
                pdbs_with_existing_loops = set([r.pdb_id for r in query])

            # get the list of pdb ids with extracted loops in data/loops
            allfiles = glob.glob(self.loop_location + "/*")
            pdbs_with_new_loops = set([f.split('/')[-1].split('_')[0] for f in allfiles])

            # skip pdb ids that already have loops:
            pdbs_needed = sorted((set(pdbs) & pdbs_with_new_loops) - pdbs_with_existing_loops)

            # process all pdb ids in pdbs that have loops in data/loops
            # once processed, the file will be moved to data/loops_done
            pdbs_needed = sorted(set(pdbs) & pdbs_with_new_loops)

            # temporary:
            # process all pdb ids that have loops in data/loops
            # once processed, the file will be moved to data/loops_done
            pdbs_needed = sorted(set(pdbs) & pdbs_with_new_loops)

        if not pdbs_needed:
            raise core.Skip("no new PDB ids need loops extracted")

        return pdbs_needed


    def query(self, session, pdb):
        """
        See if the loop_info table already has loops from this pdb
        """

        # return a query with no results when we want to process all files:
        return session.query(mod.LoopInfo).filter_by(pdb_id='XYZ')

        # skip pdb ids that already have loops:
        return session.query(mod.LoopInfo).filter_by(pdb_id=pdb)


    def remove(self, *args, **kwargs):
        """
        Does not actually remove from the database.
        We always want the loop ids to
        be consistent so we do not automatically remove loops.

        It would be good to be able to mark loops as being obsolete.
        """

        self.logger.info("We don't actually remove data from loop_info table")


    def unit_id_to_sort_name(self,loop_type,unit_ids):
        """
        Create a text string to sort the loops by
        model, symmetry, chain, residue number, and loop length
        """

        fields = loop_type.split('_')
        if fields[0][0] == "J":
            l_type = "J%02d" % int(fields[0][1:])
        else:
            l_type = fields[0] + " "

        if len(unit_ids) == 0:
            return ""

        u = unit_ids[0]
        fields = u.split('|')

        if len(fields) == 9:
            try:
                if fields[8] == "P_P":
                    sym = "P_P"
                elif "ASM" in fields[8] or "P_" in fields[8]:
                    # like ASM_2 or ASM_67 or P_23
                    sym = int(fields[8].split("_")[1])
                    sym = "%03d" % min(999,sym)
                else:
                    sym = int(fields[8].split("_")[0])
                    sym = "%03d" % min(999,sym)
            except:
                sym = '000'
        else:
            sym = '001'

        # motif length
        L = min(99,len(unit_ids))

        # loop type, model, symmetry, chain, residue number, loop length
        sort_name = "%3s%03d%2s%-4s%06d%02d" % (l_type, int(fields[1]), sym, fields[2], int(fields[4]), L)
        sort_name = sort_name.replace(" ",".")

        return sort_name

    def process_sequence(self,loop):
        """
        Loop over unit ids to pull out the sequence in a few different formats
        """
        border_total = 0

        loop["loop_name"] = ""
        loop["interior_unit_ids"] = set()
        loop["sequence_2024"] = ""
        loop["sequence_standard"] = ""
        loop["nwc_seq"] = ""

        for i, u in enumerate(loop["unit_ids"]):
            fields = u.split('|')

            a = ""
            if loop["border_indicator"][i] == "1":
                border_total += 1
                if border_total % 2 == 1:
                    # starting a new strand
                    loop["loop_name"] += "%s/%s/%s:" % (fields[1],fields[2],fields[4])
                    if i > 0 and i < len(loop["unit_ids"])-1:
                        # previous strand ended, start a new one
                        a = "*"
                else:
                    # ending a strand
                    loop["loop_name"] += fields[4]
                    if i < len(loop["unit_ids"])-1:
                        # another strand will start
                        loop["loop_name"] += ","
            else:
                loop["interior_unit_ids"].add(u)

            loop["sequence_2024"] += a
            loop["sequence_standard"] += a
            loop["nwc_seq"] += a

            # record the sequence of this nucleotide
            if fields[3] in ["A","C","G","U"]:
                s = fields[3]
                t = fields[3]
            elif fields[3] in ["DA","DC","DG","DT"]:
                s = fields[3][1]
                t = fields[3][1]
            else:
                s = "(" + fields[3] + ")"
                if fields[3] in modified_base_to_parent:
                    t = modified_base_to_parent[fields[3]]
                    if t in ['DA','DC','DG','DT']:
                        t = t[1]
                else:
                    t = s
            loop["sequence_2024"] += s
            loop["sequence_standard"] += t
            if loop["border_indicator"][i] == "0":
                loop["nwc_seq"] += s

        loop["seq"] = loop["sequence_2024"]


        if loop["loop_type"] == "IL":
            fields = loop["seq"].split("*")
            if len(fields) == 2:
                loop["r_seq"] = fields[1] + "*" + fields[0]
            else:
                print("This IL is missing a *: %s" % (str(loop)))
                self.logger.info("This IL is missing a *:\n%s" % (self.dictionary_as_text(loop)))
                return loop

            a,b = loop["nwc_seq"].split("*")
            loop["r_nwc_seq"] = b + "*" + a
        else:
            loop["r_seq"] = ""
            loop["r_nwc_seq"] = ""

        return loop


    def get_new_loops(self, pdb):
        """
        Read loops from a file created by NA_pairwise_annotations.py
        Store as a list of dictionaries with all the fields we need for
        loop_info and loop_positions tables.

        :param str pdb: PDB file to process
        :param str loop_type: The type of loop (IL, HL, J3, ...) to extract
        loops for.
        :param dict mapping: A mapping of unit ids to known loop names.
        :returns: The extracted loops.
        """

        filename = '%s_loops.txt' % pdb
        path_filename = os.path.join(self.loop_location, filename)

        """
        Format of _loops.txt file:
        J4_4TNA_001_single	4TNA|1|A|U|7,4TNA|1|A|U|8,4TNA|1|A|A|9,4TNA|1|A|2MG|10,4TNA|1|A|C|25,4TNA|1|A|M2G|26,4TNA|1|A|C|27,4TNA|1|A|G|43,4TNA|1|A|A|44,4TNA|1|A|G|45,4TNA|1|A|7MG|46,4TNA|1|A|U|47,4TNA|1|A|C|48,4TNA|1|A|5MC|49,4TNA|1|A|G|65,4TNA|1|A|A|66	1,0,0,1,1,0,1,1,0,0,0,0,0,1,1,1
        HL_4TNA_001_single	4TNA|1|A|C|13,4TNA|1|A|A|14,4TNA|1|A|G|15,4TNA|1|A|H2U|16,4TNA|1|A|H2U|17,4TNA|1|A|G|18,4TNA|1|A|G|19,4TNA|1|A|G|20,4TNA|1|A|A|21,4TNA|1|A|G|22	1,0,0,0,0,0,0,0,0,1
        """

        count = 0
        if os.path.exists(path_filename):
            with open(path_filename, 'r') as f:
                loops = f.readlines()
            loops = [l.strip() for l in loops]
            loops = [l for l in loops if l]
            count = len(loops)

        # print(path_filename)
        # print("Found %d loops in %s" % (count, pdb))

        if count == 0:
            # move the file of loops because there are none
            os.rename(os.path.join(self.loop_location, filename), os.path.join(self.done_location, filename))
            raise core.Skip("No loops found in %s" % (pdb))

        self.logger.info('Found %d loops in %s' % (count, pdb))

        # look up all basepair and stacking interactions and near basepair / stacking
        # f_lwbp_detail is not null or f_stacks is s35, s53, s33, s55 or f_lwbp is cBW or cWB
        interaction_partners = {}
        with self.session() as session:
            query = session.query(mod.UnitPairsInteractions2024.unit_id_1,mod.UnitPairsInteractions2024.unit_id_2).\
                filter(mod.UnitPairsInteractions2024.pdb_id == pdb).\
                filter(or_(mod.UnitPairsInteractions2024.f_lwbp_detail != None,
                           mod.UnitPairsInteractions2024.f_stacks != None))
            for result in query:
                if result.unit_id_1 not in interaction_partners:
                    interaction_partners[result.unit_id_1] = set()
                interaction_partners[result.unit_id_1].add(result.unit_id_2)
                if result.unit_id_2 not in interaction_partners:
                    interaction_partners[result.unit_id_2] = set()
                interaction_partners[result.unit_id_2].add(result.unit_id_1)

        # iterate over all loops in the structure
        new_loops = []
        for loop in loops:
            d = {}
            fields = loop.split('\t')
            d["loop_type"] = fields[0].split('_')[0]
            if d["loop_type"][0] == "J":
                num_strands = int(d["loop_type"][1:])
                if num_strands > 99:
                    # limit number to two digits, just in case
                    d["loop_type"] = "J00"
            d["unit_ids"]  = fields[1].split(',')
            d["unit_ids_set"]  = set(fields[1].split(','))
            d["border_indicator"] = fields[2].split(',')
            if "merged" in fields[0]:
                d["merged"] = 1
            else:
                d["merged"] = 0

            d["sort_name"] = self.unit_id_to_sort_name(d["loop_type"],d["unit_ids"])

            d = self.process_sequence(d)

            loop_unit_ids = set(d["unit_ids"])
            d["bulge"] = []
            for i, u in enumerate(d["unit_ids"]):
                fields = u.split('|')

                if not u in interaction_partners:
                    # makes no pair/stack with other nucleotides in the loop
                    d["bulge"].append(1)
                elif len(interaction_partners[u] & loop_unit_ids) == 0:
                    # makes no pair/stack with other nucleotides in the loop
                    d["bulge"].append(1)
                else:
                    d["bulge"].append(0)

            new_loops.append(d)

        # sort by model, symmetry, chain, residue number, and loop length
        new_loops = sorted(new_loops, key=lambda x: x["sort_name"])

        return new_loops, filename


    def get_existing_loops(self, pdb):
        """
        Load the loop ids that are already in the database,
        and the nucleotides in those loops.
        """

        loop_id_to_data = {}
        loop_type_to_next_id = {}

        # outerjoin was not getting lines in loop_info that were not in loop_positions
        # so just do two queries
        with self.session() as session:
            query = session.query(mod.LoopInfo.loop_id,
                                  mod.LoopInfo.loop_name,
                                  mod.LoopInfo.sort_name,
                                  mod.LoopInfo.sequence_2024,
                                  mod.LoopInfo.sequence_standard).\
                        filter(mod.LoopInfo.pdb_id == pdb)
            for result in query:
                if result.loop_id.split("_")[2] == "000":
                    continue
                if not result.loop_id in loop_id_to_data:
                    loop_id_to_data[result.loop_id] = {}
                    loop_id_to_data[result.loop_id]["unit_ids_set"] = set()
                    loop_id_to_data[result.loop_id]["interior_unit_ids"] = set()
                    loop_id_to_data[result.loop_id]["position_to_id_border"] = {}
                    loop_id_to_data[result.loop_id]["loop_name"] = result.loop_name
                    loop_id_to_data[result.loop_id]["sort_name"] = result.sort_name
                    loop_id_to_data[result.loop_id]["sequence_2024"] = result.sequence_2024
                    loop_id_to_data[result.loop_id]["sequence_standard"] = result.sequence_standard
                    loop_id_to_data[result.loop_id]["border_indicator"] = []

                fields = result.loop_id.split('_')
                loop_type = fields[0]
                loop_id_to_data[result.loop_id]["loop_type"] = loop_type

                sequential_id = int(fields[2])
                if not loop_type in loop_type_to_next_id:
                    loop_type_to_next_id[loop_type] = max(1, sequential_id+1)
                else:
                    loop_type_to_next_id[loop_type] = max(loop_type_to_next_id[loop_type], sequential_id+1)


        with self.session() as session:
            query = session.query(mod.LoopPositions.loop_id,
                                  mod.LoopPositions.position_2023,
                                  mod.LoopPositions.unit_id,
                                  mod.LoopPositions.border,
                                  mod.LoopInfo.loop_name,
                                  mod.LoopInfo.sort_name,
                                  mod.LoopInfo.sequence_2024,
                                  mod.LoopInfo.sequence_standard).\
                        outerjoin(mod.LoopInfo, mod.LoopInfo.loop_id == mod.LoopPositions.loop_id).\
                        filter(mod.LoopInfo.pdb_id == pdb)
            for result in query:
                if result.loop_id.split("_")[2] == "000":
                    continue
                if not result.loop_id in loop_id_to_data:
                    loop_id_to_data[result.loop_id] = {}
                    loop_id_to_data[result.loop_id]["unit_ids_set"] = set()
                    loop_id_to_data[result.loop_id]["interior_unit_ids"] = set()
                    loop_id_to_data[result.loop_id]["position_to_id_border"] = {}
                    loop_id_to_data[result.loop_id]["loop_name"] = result.loop_name
                    loop_id_to_data[result.loop_id]["sort_name"] = result.sort_name
                    loop_id_to_data[result.loop_id]["sequence_2024"] = result.sequence_2024
                    loop_id_to_data[result.loop_id]["sequence_standard"] = result.sequence_standard
                    loop_id_to_data[result.loop_id]["border_indicator"] = []

                loop_id_to_data[result.loop_id]["unit_ids_set"].add(result.unit_id)
                if result.border == 0:
                    loop_id_to_data[result.loop_id]["interior_unit_ids"].add(result.unit_id)
                loop_id_to_data[result.loop_id]["position_to_id_border"][result.position_2023] = (result.unit_id,result.border)

        for loop_id, old_loop in loop_id_to_data.items():
            old_loop["unit_ids"] = []
            old_loop["loop_id"] = loop_id

            border_counter = 0
            for i, data in sorted(old_loop["position_to_id_border"].items(),key = lambda x: int(x[0])):
                u,border = data
                old_loop["unit_ids"].append(u)
                old_loop["border_indicator"].append(str(border))
                border_counter += border

            # some IL with symmetry operators seem like the second strand is from model 2
            if old_loop["loop_type"] == "IL":
                if old_loop["loop_name"][0] == "1":
                    old_loop["loop_name"] = old_loop["loop_name"].replace(",2",",1")

        return loop_id_to_data, loop_type_to_next_id


    def data(self, pdb, **kwargs):
        """
        Load the loops already in the database.
        Read the loops extracted by NA_pairwise_annotations.
        Add any brand new loops
        Deprecate any old loops that should not be in the motif atlas;
        mostly ones that now have modified nucleotides in the flanking pairs.

        :param str pdb: The structure to get loops for.
        :returns: A list of all the loops.
        """

        loop_id_to_data, loop_type_to_next_id = self.get_existing_loops(pdb)

        new_loops, filename = self.get_new_loops(pdb)

        # self.logger.info("loop_type_to_next_id:")
        # self.logger.info(self.dictionary_as_text(loop_type_to_next_id))

        # accumulate data to write to tables
        loop_info = []
        loop_positions = []

        # some loop names are duplicated, and the second one won't ever get data
        loop_names_covered = set()

        for loop in new_loops:
            # default is to add the loop
            loop["add_new_loop"] = True
            exact_match_found = False

            # loop over existing loops, if any, from oldest to newest
            # I seem to have added a bunch of loops that duplicate earlier ones. Rats.
            for loop_id, old_loop in sorted(loop_id_to_data.items()):
                if not loop["loop_type"] == old_loop["loop_type"]:
                    continue

                if loop["unit_ids_set"] == old_loop["unit_ids_set"] or loop["loop_name"] == old_loop["loop_name"]:
                    # exact match
                    exact_match_found = True
                    loop_names_covered.add(old_loop["loop_name"])
                    loop_names_covered.add(loop["loop_name"])

                    if loop["unit_ids_set"] == old_loop["unit_ids_set"]:
                        self.logger.info('Exact unit_id match between old %s with new %s' % (loop_id, loop["loop_name"]))
                    else:
                        self.logger.info('Exact loop_name match between old %s with new %s' % (loop_id, loop["loop_name"]))

                    # do not add the new loop itself, instead update the old loop where necessary
                    loop["add_new_loop"] = False

                    if old_loop["sort_name"] == "" or old_loop["sequence_2024"] == None or not old_loop["sequence_2024"] == loop["sequence_2024"]:
                        # old loop was not processed on a previous run, or not processed correctly
                        # check and update the information for the loop_positions table
                        for i, u in enumerate(loop["unit_ids"]):
                            if not i+1 in old_loop["position_to_id_border"]:
                                # You can only add when the position is not already in the table
                                # Just skip these, even though now some unit ids may be wrong
                                # There is only so much you can do with old messed up loops
                                # or old_loop["position_to_id_border"][i][0] != u \
                                # or old_loop["position_to_id_border"][i][1] != loop["border_indicator"][i]:

                                self.logger.info('Backfilled loop_positions: %s %3d %s' % (loop_id, i+1, u))
                                loop_positions.append(mod.LoopPositions(
                                    loop_id = loop_id,
                                    position = i+1,
                                    bulge = loop["bulge"][i],
                                    flanking = loop["border_indicator"][i],
                                    border = loop["border_indicator"][i],
                                    unit_id = u,
                                    position_2023 = i+1))

                        # many old loops have incorrect length, seq, r_seq, nwc_seq, r_nwc_seq
                        # because modified nucleotides were ignored.
                        # fix them by copying the new information to the old loop
                        loop_id_to_data[loop_id] = loop
                        loop_id_to_data[loop_id]["copied"] = True
                        if not "deprecate" in loop_id_to_data[loop_id]:
                            loop_id_to_data[loop_id]["deprecate"] = 0

                    else:
                        # we have already updated this loop on a previous pass
                        old_loop["save_loop_info"] = False

                    break
                elif loop["interior_unit_ids"] & old_loop["interior_unit_ids"]:
                    print('Partial match between old %s and new %s %s' % (loop_id, loop["loop_type"], loop["loop_name"]))
                    self.logger.info('Partial match between old %s with %s' % (loop_id, loop["loop_name"]))
                    self.logger.info('  https://rna.bgsu.edu/rna3dhub/loops/view/%s' % loop_id)
                    if loop["merged"] == 1:
                        print('  New loop is merged')
                        self.logger.info('  New loop is merged, so add it')
                        self.logger.info('  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s' % ",".join(loop["unit_ids"]))

                        if not "deprecate" in old_loop:
                            old_loop["deprecate"] = 0

                    elif loop["unit_ids_set"] <= old_loop["unit_ids_set"]:
                        print("  New loop is entirely contained in the old loop")
                        print("  Old loop also has these nucleotides: %s" % (old_loop["unit_ids_set"] - loop["unit_ids_set"]))
                        self.logger.info("  New loop is entirely contained in the old loop")
                        self.logger.info("  Old loop also has these nucleotides: %s" % (old_loop["unit_ids_set"] - loop["unit_ids_set"]))
                        self.logger.info("  Deprecate the old loop")

                        old_loop["deprecate"] = 1

                    elif old_loop["unit_ids_set"] <= loop["unit_ids_set"]:
                        print("  Old loop is entirely contained in the new loop")
                        print("  New loop also has these nucleotides: %s" % (loop["unit_ids_set"] - old_loop["unit_ids_set"]))
                        self.logger.info("  Old loop is entirely contained in the new loop")
                        self.logger.info("  New loop also has these nucleotides: %s" % (loop["unit_ids_set"] - old_loop["unit_ids_set"]))
                        self.logger.info("  Bigger loop now, why?")
                    else:
                        print("  Old loop also contains: %s" % (old_loop["unit_ids_set"] - loop["unit_ids_set"]))
                        print("  New loop also contains: %s" % (loop["unit_ids_set"] - old_loop["unit_ids_set"]))
                        self.logger.info("  Complicated overlap")
                        self.logger.info("  Old loop also contains: %s" % (old_loop["unit_ids_set"] - loop["unit_ids_set"]))
                        self.logger.info("  New loop also contains: %s" % (loop["unit_ids_set"] - old_loop["unit_ids_set"]))

            if exact_match_found:
                # done processing this new loop, move on to the next new loop
                continue

        # accumulate database entries for each new loop
        for loop in new_loops:
            if loop.get("add_new_loop",False):

                loop_names_covered.add(loop["loop_name"])

                # self.logger.info("loop_type_to_next_id:")
                # self.logger.info(self.dictionary_as_text(loop_type_to_next_id))

                loop_type = loop["loop_type"]
                if loop_type not in loop_type_to_next_id:
                    loop_type_to_next_id[loop_type] = 1
                id = loop_type_to_next_id[loop_type]
                if id <= 999:
                    sequential_id = '%03d' % id
                else:
                    sequential_id = '%06d' % id
                loop_id = '%s_%s_%s' % (loop_type, pdb, sequential_id)
                loop_type_to_next_id[loop_type] += 1

                loop_info.append(mod.LoopInfo(
                    loop_id=loop_id,
                    type=loop_type,
                    pdb_id=pdb,
                    sequential_id=sequential_id,
                    length=len(loop["unit_ids"]),
                    seq=loop["seq"],
                    r_seq=loop["r_seq"],
                    nwc_seq=loop["nwc_seq"],
                    r_nwc_seq=loop["r_nwc_seq"],
                    loop_name=loop["loop_name"],
                    sequence_2024=loop["sequence_2024"],
                    sequence_standard=loop["sequence_standard"],
                    sort_name=loop["sort_name"],
                    merged=loop["merged"]))

                for i, u in enumerate(loop["unit_ids"]):
                    self.logger.info('Brand new loop_positions: %s %s %s' % (loop_id, i+1, u))
                    loop_positions.append(mod.LoopPositions(
                        loop_id=loop_id,
                        position = i+1,
                        bulge = loop["bulge"][i],
                        flanking = loop["border_indicator"][i],
                        border = loop["border_indicator"][i],
                        unit_id = u,
                        position_2023 = i+1))

                loop["done"] = True

        self.logger.info('Adding %d brand new entries to loop_info' % len(loop_info))
        for d in loop_info:
            self.logger.info(self.dictionary_as_text(d.__dict__))

        # update information about old loops for loop_info table
        # presume that loop_positions is correct when an exact match is not available
        old_loop_info = []
        for loop_id, old_loop in list(loop_id_to_data.items()):
            if old_loop.get("save_loop_info",True) == False:
                # this loop was already processed, do not save data again
                continue

            if old_loop.get("copied",False):
                # data was already fixed by copying the exactly matching new loop
                pass
            else:
                # these are old loops that are not exact matches
                # some of them are missing loop_position data
                # some of them should be deleted because they are not really loops,
                # but that is hard to determine automatically, so don't worry about them

                old_loop["sort_name"] = self.unit_id_to_sort_name(loop_id,old_loop["unit_ids"])

                if len(old_loop["sort_name"]) == 0:
                    # no unit ids, nothing useful to do here
                    self.logger.info("  No unit ids for %s" % loop_id)
                    continue

                # the next line may fail if loop_positions is missing data or is wrong
                try:
                    old_loop = self.process_sequence(old_loop)
                except:
                    self.logger.info("  Unable to process sequence for %s" % loop_id)
                    self.logger.info("  select * from loop_positions where loop_id = '%s';" % loop_id)
                    self.logger.info(self.dictionary_as_text(old_loop))
                    continue

                if old_loop["loop_name"] == "":
                    self.logger.info("  No loop name for %s" % loop_id)
                    continue

                # deprecate some problematic "loops" that are not really loops
                if old_loop["loop_type"] == "HL" and len(old_loop["unit_ids"]) == 2:
                    old_loop["deprecate"] = 1
                    print("Deprecating %s with two nucleotides" % loop_id)
                    self.logger.info("Deprecating %s with two nucleotides" % loop_id)
                elif old_loop["loop_type"] == "IL" and len(old_loop["unit_ids"]) == 4:
                    old_loop["deprecate"] = 1
                    print("Deprecating %s with four nucleotides" % loop_id)
                    self.logger.info("Deprecating %s with four nucleotides" % loop_id)
                elif old_loop["loop_type"] == "J3" and len(old_loop["unit_ids"]) == 6:
                    old_loop["deprecate"] = 1
                    print("Deprecating %s with six nucleotides" % loop_id)
                    self.logger.info("Deprecating %s with six nucleotides" % loop_id)

            loop_names_covered.add(old_loop["loop_name"])

            old_loop_info.append(mod.LoopInfo(
                loop_id=loop_id,
                length=len(old_loop["unit_ids"]),
                r_seq=old_loop["r_seq"],
                nwc_seq=old_loop["nwc_seq"],
                r_nwc_seq=old_loop["r_nwc_seq"],
                seq=old_loop["seq"],
                sequence_2024=old_loop["sequence_2024"],
                sequence_standard=old_loop["sequence_standard"],
                sort_name=old_loop["sort_name"],
                deprecate=old_loop.get("deprecate",0)))

        self.logger.info('Adding %d updated entries to loop_info from old loops that were not exact matches' % len(old_loop_info))
        for d in old_loop_info:
            self.logger.info('Updated loop_info entry:')
            self.logger.info(self.dictionary_as_text(d.__dict__))


        # move the file of loops once all loops have been processed
        all_loops_done = True
        for loop in new_loops:
            if not loop.get("done",False) and not loop.get("loop_name","no") in loop_names_covered:
                all_loops_done = False
                break
        for loop_id, old_loop in loop_id_to_data.items():
            if not old_loop.get("done",False) and not old_loop.get("loop_name","no") in loop_names_covered:
                all_loops_done = False
                break
        if all_loops_done:
            self.logger.info('All loops processed, moving %s to %s' % (filename, self.done_location))
            # comment out the next two lines for testing
            os.rename(os.path.join(self.loop_location, filename), os.path.join(self.done_location, filename))
        return loop_info + loop_positions + old_loop_info

        # while testing, do not actually add the new loops to the database
        print('Not writing to the database')
        return []
