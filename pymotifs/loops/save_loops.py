"""
Read HL, IL, J3, etc. from a file created by NA_pairwise_annotations
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
    loop_types = ['IL', 'HL', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'J11', 'J12', 'J13', 'J14', 'J15']
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
            print(key, d[key])


    def to_process(self, pdbs, **kwargs):

        # get the list of pdb ids that already have loops in the database
        with self.session() as session:
            query = session.query(mod.LoopInfo.pdb_id).distinct()
            pdbs_with_existing_loops = set([r.pdb_id for r in query])

        # get the list of pdb ids with extracted loops in data/loops
        allfiles = glob.glob(self.loop_location + "/*")
        pdbs_with_new_loops = set([f.split('/')[-1].split('_')[0] for f in allfiles])

        pdbs_needed = sorted((set(pdbs) & pdbs_with_new_loops) - pdbs_with_existing_loops)


        if False:
            # code to delete newly made loops but not overlap with the motif atlas

            # query the database to find all loop ids in loop_info where the length of
            # the field sort_name is greater than 0
            with self.session() as session:
                query = session.query(mod.LoopInfo.loop_id).\
                    filter(mod.LoopInfo.sort_name != '')

                saved_loop_ids = set([r.loop_id for r in query])

            print('Found %5d loops in loop_info with non-trivial sort_name' % len(saved_loop_ids))

            # query the motif atlas for loop ids that are in the motif atlas
            with self.session() as session:
                query = session.query(mod.MlLoops.loop_id).distinct()
                motif_atlas_loop_ids = set([r.loop_id for r in query])

            print('Found %5d loops in motif_atlas' % len(motif_atlas_loop_ids))
            print('Found %5d loops in both' % len(saved_loop_ids & motif_atlas_loop_ids))
            print('They are: %s' % (saved_loop_ids & motif_atlas_loop_ids))


        if not pdbs_needed:
            raise core.Skip("no new PDB ids need loops extracted")

        return pdbs_needed


    def query(self, session, pdb):
        """
        See if the loop_info table already has loops from this pdb
        """

        # return a query with no results when we want to process all:
        # return session.query(mod.LoopInfo).filter_by(pdb_id='XYZ')

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


    def _next_loop_number_string(self, current):
        """
        Compute the next loop number string. This will pad to either 3 or 6
        characters with zeros. If the next number is over 999 we use 6,
        otherwise 3 as the length to pad to.

        :current: The current loop count.
        :returns: A string of the next loop id.
        """

        next_number = current + 1
        if next_number > 999:
            return str(next_number).rjust(6, '0')
        return str(next_number).rjust(3, '0')


    def _get_loop_id(self, units, pdb_id, loop_type, mapping):
        """
        Compute the loop id to use for the given unit string. This will
        build a string like IL_1S72_001 or IL_4V4Q_001.
        In structures with over 999 loops, we will pad with zeros to 6
        characters, but keep the standard padding to 3 characters otherwise.

        :units: The concatenated unit or nt id string.
        :pdb_id: The pdb id to use.
        :loop_type: The type of loop.
        :mapping: A mapping from unit string to known loop_id.
        :returns: A string of the new loop id.
        """

        if units not in mapping:
            str_number = self._next_loop_number_string(len(mapping))
            loop_id = '%s_%s_%s' % (loop_type, pdb_id, str_number)
            self.logger.info('Created new loop id %s, for nucleotides %s',
                             loop_id, units)
            mapping[units] = str(loop_id)

        return mapping[units]


    # def _get_fake_loop_id(self, pdb_id, loop_type):
    #     """
    #     Indicate that a structure has been checked for a certain type
    #     of loop and none were found.

    #     We won't need this for Python-extracted loops, because
    #     we will extract loops when we do pairwise annotations, so we
    #     won't ever try to extract loops from a structure that has none.
    #     """

    #     loop_id = '%s_%s_%s' % (loop_type, pdb_id, '000')
    #     self.logger.info('Created new fake loop id %s for pdb_id %s', loop_id, pdb_id)

    #     return loop_id


    def get_existing_loops(self, pdb):
        """
        Load the loop ids that are already in the database,
        and the nucleotides in those loops.
        """

        old_loops = []
        loop_type_to_next_id = {}
        for loop_type in self.loop_types:
            loop_type_to_next_id[loop_type] = 1

        # with self.session() as session:
        #     query = session.query(mod.UnitPairsInteractions2024.unit_id_1,mod.UnitPairsInteractions2024.unit_id_2).\
        #         filter(mod.UnitPairsInteractions2024.pdb_id == pdb).\
        #         filter(or_(mod.UnitPairsInteractions2024.f_lwbp_detail != None,
        #                    mod.UnitPairsInteractions2024.f_stacks != None))
        #     for result in query:


        return old_loops, loop_type_to_next_id


    def get_new_loops(self, pdb):
        """
        Read loops from a file created by NA_pairwise_annotations
        Store as a list of dictionaries with all the fields we need for
        loop_info and loop_positions tables.

        Uses matlab to extract the loops for a given structure of a specific
        type. This will also save the loop files into the correct place.

        :param str pdb: PDB file to process
        :param str loop_type: The type of loop (IL, HL, J3, ...) to extract
        loops for.
        :param dict mapping: A mapping of unit ids to known loop names.
        :returns: The extracted loops.
        """

        filename = '%s_loops.txt' % pdb
        path_filename = os.path.join(self.loop_location, filename)

        """
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
            d["unit_ids"]  = fields[1].split(',')
            d["border_indicator"] = fields[2].split(',')
            d["sequence_2024"] = ""
            d["sequence_standard"] = ""
            d["nwc_seq"] = ""
            d["loop_name"] = ""
            d["bulge"] = []
            if "merged" in fields[0]:
                d["merged"] = 1
            else:
                d["merged"] = 0

            border_total = 0
            loop_unit_ids = set(d["unit_ids"])

            for i, u in enumerate(d["unit_ids"]):
                fields = u.split('|')

                if len(fields) == 9:
                    sym = fields[8][0]
                else:
                    sym = '1'

                if i == 0:
                    # enable sorting entire loops by model, symmetry, chain, first residue number, and loop length
                    d["sort_name"] = "%03d_%s_%-4s_%06d_%03d" % (int(fields[1]), sym, fields[2], int(fields[4]), len(d["unit_ids"]))
                    d["sort_name"] = d["sort_name"].replace(" ",".")

                a = ""
                if d["border_indicator"][i] == "1":
                    border_total += 1
                    if border_total % 2 == 1:
                        # starting a new strand
                        d["loop_name"] += "%s/%s/%s:" % (fields[1],fields[2],fields[4])
                        if i > 0 and i < len(d["unit_ids"])-1:
                            # previous strand ended, start a new one
                            a = "*"
                    else:
                        # ending a strand
                        d["loop_name"] += fields[4]
                        if i < len(d["unit_ids"])-1:
                            # another strand will start
                            d["loop_name"] += ","

                d["sequence_2024"] += a
                d["sequence_standard"] += a
                d["nwc_seq"] += a

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
                d["sequence_2024"] += s
                d["sequence_standard"] += t
                if d["border_indicator"][i] == "0":
                    d["nwc_seq"] += s

                if not u in interaction_partners:
                    # makes no pair/stack with other nucleotides in the loop
                    d["bulge"].append(1)
                elif len(interaction_partners[u] & loop_unit_ids) == 0:
                    # makes no pair/stack with other nucleotides in the loop
                    d["bulge"].append(1)
                else:
                    d["bulge"].append(0)

            d["seq"] = d["sequence_2024"]

            # self.print_dictionary(d)
            # print("")

            if d["loop_type"] == "IL":
                a,b = d["seq"].split("*")
                d["r_seq"] = b + "*" + a
                a,b = d["nwc_seq"].split("*")
                d["r_nwc_seq"] = b + "*" + a
            else:
                d["r_seq"] = ""
                d["r_nwc_seq"] = ""

            new_loops.append(d)

        # sort by model, chain, residue number, and loop length
        # this will put merged loops in their correct position in the structure
        new_loops = sorted(new_loops, key=lambda x: x["sort_name"])

        return new_loops, filename


    def _mapping(self, pdb_id, loop_type, normalizer):
        """
        Compute a mapping from the nts to the loop id.
        This is used for setting ids by either looking up the old
        known id or creating a new one if no one is found.

        :pdb_id: The pdb id to search.
        :loop_type: The loop type.
        :returns: A dictonary mapping from the nts to the loop id. If there are
        no loops it returns an empty dictonary.
        """

        mapping = {}
        with self.session() as session:
            query = self.query(session, pdb_id).filter_by(type=loop_type)
            for result in query:
                unit_ids = normalizer(result.unit_ids)
                if unit_ids in mapping:
                    self.logger.error("Loop %s duplicates %s",
                                      result.loop_id, mapping[unit_ids])
                    continue
                mapping[unit_ids] = result.loop_id
        return mapping


    def data(self, pdb, **kwargs):
        """
        Load the loops already in the database.
        Read the loops extracted by NA_pairwise_annotations.
        Add any brand new loops
        Deprecate any old loops that should be replaced with new loops ... in the future

        :param str pdb: The structure to get loops for.
        :returns: A list of all the loops.
        """

        old_loops, loop_type_to_next_id = self.get_existing_loops(pdb)

        # # only process new structures
        # if len(old_loops) > 0:
        #     raise core.Skip("Loops already exist for %s" % pdb)

        new_loops, filename = self.get_new_loops(pdb)

        deprecate_loops = []
        partial_matches = []
        for loop in new_loops:
            # default is to add the loop
            loop["add_loop"] = True

            exact_match_found = False
            for old_loop in old_loops:
                if loop["unit_ids_set"] == old_loop["unit_ids_set"]:
                    exact_match_found = True
                elif loop["unit_ids_set"] & old_loop["unit_ids_set"]:
                    partial_matches.append((loop, old_loop))

            if exact_match_found:
                # nothing to do
                loop["add_loop"] = False
                continue

            if len(partial_matches) > 0:
                # deprecate the old loop
                # make a new loop with a brand new loop id
                pass


        # will need to add sort_name to every loop so we can sort by that

        # accumulate database entries for each new loop
        loop_info = []
        loop_positions = []
        for loop in new_loops:
            if loop["add_loop"]:
                loop_type = loop["loop_type"]
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
                    loop_positions.append(mod.LoopPositions(
                        loop_id=loop_id,
                        position = i+1,
                        bulge = loop["bulge"][i],
                        flanking = loop["border_indicator"][i],
                        border = loop["border_indicator"][i],
                        unit_id = u,
                        position_2023 = i))

        self.logger.info('Adding %d new entries to loop_info' % len(new_loops))
        self.logger.info('Adding %d new entries to loop_positions' % len(loop_positions))


        # note that the file of loops was already processed
        os.rename(os.path.join(self.loop_location, filename), os.path.join(self.done_location, filename))

        return loop_info + loop_positions
