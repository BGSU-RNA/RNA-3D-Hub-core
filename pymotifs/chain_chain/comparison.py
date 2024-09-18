"""
Compute chain to chain discrepancies within groups having sequence correspondences.
This will:
1) look at good correspondences,
2) extract all the aligned chains,
3) compute the geometric discrepancy between them, and
4) place the discrepanices
in the database.

This will only compare the first chain in each IFE.
"""

import functools as ft
import itertools as it
import numpy as np
import operator as op
import os
import pickle

from collections import defaultdict

from sqlalchemy import and_
from sqlalchemy import or_
from sqlalchemy.orm import aliased
# from sqlalchemy.sql import union_all

from fr3d.geometry.discrepancy import matrix_discrepancy

from pymotifs import core
import pymotifs.utils as ut
from pymotifs import models as mod
from pymotifs.utils import discrepancy as disc

from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqUnitMappingLoader
from pymotifs.ife.loader import Loader as IfeLoader
# from pymotifs.units.centers import Loader as CenterLoader
# from pymotifs.units.rotation import Loader as RotationLoader
from pymotifs.units.center_rotation import Loader as CenterRotationsLoader

from pymotifs.nr.groups.simplified import Grouper

def pick(preferences, key, iterable):
    """Pick the most preferred value from a list of possibilities.

    Parameters
    ----------
    preferences : list
        A list of possibilities to select from.
    key : str
        Key to define the attribute
    iterable : iterable
        The iterable to get the unique values from.

    Returns
    -------
    choice : obj
        A member of preferences if any exist, otherwise the alphabetically
        first choice from the iterable.
    """

    possible = set(getattr(result, key) for result in iterable)
    if not possible:
        raise core.InvalidState("Nothing to pick from")

    for pref in preferences:
        if pref in possible:
            return pref
    return sorted(possible)[0]


class Loader(core.SimpleLoader):
    """A Loader to get all chain to chain similarity data. This will use the
    correspondences between ifes to compute the discrepancy between them.
    """

    """The tuple of chains can't be marked in the database."""
    mark = False

    """We allow for no data to be written when appropriate."""
    allow_no_data = True

    """The dependencies for this stage"""
    dependencies = set([CorrespondenceLoader, ExpSeqUnitMappingLoader,
                        IfeLoader, CenterRotationsLoader])
    # dependencies = set([])


    def known_unit_entries(self, table, molecule_type):
        """
        Create a set of (pdb, chain) tuples for all chains that have entries
        in the given table. This relies upon the table having a unit_id column
        that can be joined to unit_info.

        Parameters
        ----------
        table : Table
            The table to join against.

        Returns
        -------
        known : set
            A set of tuples of (pdb, chain).
        """

        # filtering by unit_type_id might make this much slower ... paradoxically

        with self.session() as session:
            info = mod.UnitInfo
            query = session.query(info.pdb_id,
                                  info.chain,
                                  ).\
                join(table, table.unit_id == info.unit_id).\
                filter(info.unit_type_id  == molecule_type).\
                distinct()
            return set((r.pdb_id, r.chain) for r in query)


    def is_member(self, known, chain):
        """Check if a chain is a member of the set of knowns.

        Parameters
        ----------
        known : set
            A set of tuples as produced by known_unit_entries.
        chain : dict
            A dict with 'pdb' and 'name' entries for the pdb id and chain name.

        Returns
        -------
        member : bool
            True if the chain is a member of the set.
        """
        getter = op.itemgetter('pdb', 'name')
        return getter(chain) in known


    def to_process(self, pdbs, **kwargs):
        """
        This will compute all pairs to compare. This will group all pdbs
        using only sequence and species and then produce a list of chains that
        are only the chains in the same group. It will also filter out all
        pairs that have a chain with a poor resolution or is too small. This
        will produce all comparisons that are needed for the NR set and no
        more. If no groups are produced we raise an exception.

        Parameters
        ----------
        pdbs : list
            List of pdb ids to group by species and sequence

        Raises
        ------
        pymotifs.core.InvalidState
            If there is no grouping produced.

        Returns
        -------
        pairs : list
            A list of (first, seconds) pairs of chain ids to compare.
        """

        if len(pdbs) < 100:
            self.logger.info('Not running chain_chain.comparison for small number of pdbs')
            raise core.Skip("Not enough pdbs to compare")

        pdbs_set = set(pdbs)

        # while developing DNA equivalence classes, run those manually
        nr_molecule_parent_current = kwargs.get('nr_molecule_parent_current','')

        if nr_molecule_parent_current and 'dna' in nr_molecule_parent_current.lower():
            molecule_types = ['dna']
        else:
            molecule_types = ['rna']

        # create lists of chain ids from each group, then data method will loop over all groups of chain ids
        groups_of_chain_ids = []

        for molecule_type in molecule_types:

            self.logger.info('Processing molecule type %s' % molecule_type)

            # Group PDB ids by species and sequence
            # groups is a list of lists of integer chain identifiers
            # grouper is part of the nr stage, so it appears that way in the log file
            grouper = Grouper(self.config, self.session)
            grouper.use_discrepancy = False                                                                                      ## change the use_discrepancy to False.
            grouper.must_enforce_single_species = False   # use False with DNA
            grouper.must_enforce_single_species = True    # use True on production
            grouper.molecule_type = molecule_type
            groups = grouper(pdbs)                                                                                               ## Grouped PDB ids by sequence and species. How this class works
            if not groups:
                raise core.InvalidState("No groups produced")

            self.logger.info('Found %d groups of IFEs of type %s' % (len(groups),molecule_type))

            found_non_singleton_group = False
            for group in groups:
                if 'members' in group and len(group['members']) > 1:
                    found_non_singleton_group = True
                    break
            if not found_non_singleton_group:
                raise core.Skip("No groups of multiple chains to compare")

            self.logger.info('Setting up has_rotations for the given molecule type')

            # only keep chain ids where the chain has at least one nucleotide with a
            # base center and rotation matrix
            # or we could just wait and look that up later
            # this query is slow, it would be nice to avoid it
            has_rotations = ft.partial(self.is_member,
                                    self.known_unit_entries(mod.UnitRotations,molecule_type))                                                  ## modify `is_member` function

            # it's enough to check for rotations; if no rotations, then no centers
            # has_centers = ft.partial(self.is_member,
            #                         self.known_unit_entries(mod.UnitCenters))

            # the variable group['members'] is a list of dictionaries with tons of information

            self.logger.info("Filtering out chains without centers and rotations")
            for group in groups:
                self.logger.info('Checking a group of size %d' % len(group['members']))
                chains = filter(disc.valid_chain, group['members'])                                                           ##
                chains = filter(has_rotations, chains)
                # chains = filter(has_centers, chains)
                chains = map(op.itemgetter('db_id'), chains)
                # convert chains from iterator to list of chain ids and append to list
                # use a set to not repeat any chain ids, was a problem with 2M4Q|1|1
                chain_list = sorted(set(chains))
                if len(chain_list) > 1:
                    groups_of_chain_ids.append(chain_list)

                if len(chain_list) < len(group['members']):
                    self.logger.info('Dropped from %d to %d members in group' % (len(group['members']),len(chain_list)))


        # Note how many groups are left now that we filtered as above
        self.logger.info("Found %d groups with at least two chains" % len(groups_of_chain_ids))

        # start with the smallest group
        # groups_of_chain_ids.sort(key=len)

        # start with the largest group
        groups_of_chain_ids.sort(key=len,reverse=True)

        if len(groups_of_chain_ids) == 0:
            raise core.Skip("No groups of chains to compare")

        GeneratePickleFiles = False   # appropriate to use when debugging the rest of the program
        GeneratePickleFiles = True    # must be used in production, to update the files each week

        if len(pdbs) < 500:
            GeneratePickleFiles = False
            self.logger.info('Not generating pickle files for small number of pdbs')

        if GeneratePickleFiles:
            self.logger.info('Getting the unit_id to experimental sequence position table')
            # query and write to disk unit correspondence data once per run of chain_chain/comparison.py
            # takes about 2 1/3 minutes on production in December 2020
            # takes over 15 minutes when re-written to check sym_op, alt_id, glycosidic center

            known_pdbs = set()
            chain_unit_to_position = defaultdict(dict)

            # read previous chain unit to position data, takes 45 seconds on production on 2024-06-17
            if os.path.exists("unit_to_position.pickle"):
                self.logger.info('Reading pickle file unit_to_position.pickle')
                chain_unit_to_position = pickle.load(open("unit_to_position.pickle", "rb" ))
                for chain in chain_unit_to_position.keys():
                    known_pdbs.add(chain.split("|")[0])
                self.logger.info('Read mappings for %d pdb files from unit_to_position.pickle' % len(known_pdbs))

            needed_pdbs = set(pdbs) - known_pdbs

            # the query below is slow when being run on 8000 pdbs, which would happen
            # if the .pickle file above is deleted or corrupted
            # mysql show processlist indicates that mysql uses a temporary table on disk
            # on 2024-07-24 it took 22 minutes
            # the query below could be split to do 1000 pdbs at a time

            self.logger.info('Updating unit_to_position.pickle for %d pdbs' % len(needed_pdbs))
            self.logger.info(sorted(needed_pdbs))

            # get correspondences between unit ids and experimental sequence positions
            # sort to put no alt_id before A before B before C, etc.
            with self.session() as session:
                EM = mod.ExpSeqUnitMapping
                UI = mod.UnitInfo
                UC = mod.UnitCenters
                query = session.query(EM.unit_id,EM.exp_seq_position_id,UI.sym_op).\
                    join(UI, UI.unit_id == EM.unit_id).\
                    join(UC, UC.unit_id == EM.unit_id).\
                    filter(UC.name == 'glycosidic').\
                    filter(UI.pdb_id.in_(needed_pdbs)).\
                    order_by(UI.pdb_id,UI.model,UI.chain,UI.alt_id.is_(None).desc(),UI.alt_id)

                chain_to_symmetries = defaultdict(set)
                count = 0
                for r in query:
                    # self.logger.info('Unit id %s, position %s, sym_op %s' % (r.unit_id,r.exp_seq_position_id,r.sym_op))
                    count += 1
                    if r.unit_id and "|" in r.unit_id:    # not sure why, but some rows have None
                        fields = r.unit_id.split("|")
                        if len(fields) > 3:
                            ## skip unexpected pdbs
                            if not fields[0] in needed_pdbs:
                                continue
                            chain = "|".join(fields[0:3])   # pdb id, model, chain is the top level key
                            chain_to_symmetries[chain].add(r.sym_op)
                self.logger.info('Got %d raw unit to position mappings' % count)

                for chain, symmetries in chain_to_symmetries.items():
                    if '1_555' in symmetries:
                        # use the default symmetry when available
                        symmetry = '1_555'
                    elif len(symmetries) >= 1:
                        # otherwise use the lowest numbered one, whatever that means
                        symmetry = sorted(symmetries)[0]
                    else:
                        symmetry = ''
                        self.logger.info('No symmetry operators for %s' % chain)
                    # record the symmetry to use for this chain
                    chain_to_symmetries[chain] = symmetry

                chain_to_simple_unit_id = defaultdict(set)
                count = 0
                # looping over the query again may be just as slow as the first time
                for r in query:
                    if r.unit_id and "|" in r.unit_id:    # not sure why, but some rows have None
                        fields = r.unit_id.split("|")
                        if len(fields) >= 5:
                            ## skip unexpected pdbs, just in case
                            if not fields[0] in pdbs_set:
                                continue

                            chain = "|".join(fields[0:3])   # pdb id, model, chain is the top level key
                            # only store unit ids from the one designated symmetry operator
                            if r.sym_op == chain_to_symmetries[chain]:
                                lf = len(fields)
                                if lf in [5,6,7]:
                                    # plain or just an alt id
                                    simple_unit_id = "|".join([fields[0],fields[1],fields[2],'',fields[4]])  # remove sequence
                                elif lf == 8:
                                    # insertion code, possibly with alt_id
                                    # remove sequence and alt_id
                                    simple_unit_id = "|".join([fields[0],fields[1],fields[2],'',fields[4],fields[5],'',fields[7]])
                                else:
                                    # symmetry, possibly with insertion code, possibly with alt_id
                                    # remove sequence and alt_id
                                    simple_unit_id = "|".join([fields[0],fields[1],fields[2],'',fields[4],fields[5],'',fields[7],fields[8]])

                                # only store one version of each unit id
                                if not simple_unit_id in chain_to_simple_unit_id[chain]:
                                    chain_to_simple_unit_id[chain].add(simple_unit_id)
                                    chain_unit_to_position[chain][r.unit_id] = r.exp_seq_position_id
                                    count += 1
                                else:
                                    self.logger.info('Unit id %s would duplicate earlier %s' % (r.unit_id,simple_unit_id))

                self.logger.info('Added %d final unit to position mappings' % count)
                pickle.dump(chain_unit_to_position, open("unit_to_position.pickle", "wb" ), 2)
                self.logger.info('Wrote pickle file unit_to_position.pickle')

            del chain_unit_to_position     # clear this from memory

        if GeneratePickleFiles:
            # took about 6.5 minutes on production in December 2020
            # took about 22 minutes on production in June 2024
            # Maybe it would be faster to update incrementally, only getting new correspondences
            self.logger.info('Getting the position to position mappings')

            position_to_position = defaultdict(list)                                                            ## A dictionary-like object. It will not raise error. It will provides a default value for the key does not exist. In this case, it will return a empty list object if the key does not exist.
            count = 0

            i = 0
            j = 1
            chunk_size = 100000
            newcount = 1
            found_lines = False

            while j < 9999:
                if newcount == 0 and found_lines:
                    # already found data, but then got an empty set of data
                    j = 9999            # this will be the last query, check for very large numbers
                else:
                    j = i+1

                with self.session() as session:
                    CP = mod.CorrespondencePositions
                    query = session.query(CP.exp_seq_position_id_1,CP.exp_seq_position_id_2).\
                    select_from(CP).\
                    filter(CP.exp_seq_position_id_1 < chunk_size*j).\
                    filter(CP.exp_seq_position_id_1 >= chunk_size*i)

                    newcount = 0
                    for r in query:
                        position_to_position[r.exp_seq_position_id_1].append(r.exp_seq_position_id_2)                       ## the result is like: defaultdict(list, {'aaa': ['bbb', 'fff'], 'ccc': ['ddd']})
                        count += 1
                        newcount += 1

                    print('Total of %10d position to position correspondences, id up to %10d' % (count,chunk_size*j))
                    self.logger.info('Total of %10d position to position correspondences, id up to %10d' % (count,chunk_size*j))

                    if newcount > 0:
                        found_lines = True

                i = i + 1

            pickle.dump(position_to_position, open("position_to_position.pickle", "wb" ), 2)                                ## save the position_to_position object
            self.logger.info('Wrote pickle file position_to_position.pickle')
            self.logger.info('Maximum exp_seq_position_id_1 is %d' % max(position_to_position.keys()))

            del position_to_position   # clear the memory

        return groups_of_chain_ids


    def is_missing(self, entry, **kwargs):
        """Determine if we do not have any data. If we have no data then we
        will recompute. This method must be implemented by inheriting classes
        and is how we determine if we have data or not.

        Parameters
        ----------
        entry : object
            The data to check
        **kwargs : dict
            Generic keyword arguments

        Returns
        -------
        missing : bool
            True if the data is missing.
        """
        return True


    def query(self, session, pair):
        """Create a query to find the comparisons using the given pair of
        chains. This will find a pair in either direction. Also, this ignores
        the model and other information and uses only chain ids.

        Parameters
        ----------
        session : pymotifs.core.Session
            The session to use

        pair : (int, int)
            The pair of chain ids to look up.

        Returns
        -------
        query : query
            The query for chains.
        """

        self.logger.info("query: first: %s // second (clean): %s" % (pair[0], pair[1][0]))

        sim = mod.ChainChainSimilarity

        return session.query(sim).\
            filter(or_(and_(sim.chain_id_1==pair[0],sim.chain_id_2.in_(pair[1])),
                       and_(sim.chain_id_2==pair[0],sim.chain_id_1.in_(pair[1]))))

    def get_unit_correspondences_intersect(self,info1,info2,unit_to_position,position_to_position):
        """
        Given two chains, use unit_to_position and position_to_position
        mappings to find unit to unit correspondences.
        Return a list of pairs of units
        """

        matching_pairs = []

        # for now, concentrate on the first chain in any multi-chain IFEs
        chain1 = info1['ife_id'].split('+')[0]
        chain2 = info2['ife_id'].split('+')[0]

        # number of resolved nucleotides
        length1 = len(unit_to_position.get(chain1,[]))
        length2 = len(unit_to_position.get(chain2,[]))

        if length1 == 0:
            self.logger.warning("No unit to position mapping for %s" % chain1)
        elif length2 == 0:
            self.logger.warning("No unit to position mapping for %s" % chain2)
        elif length1/length2 > 10 or length2/length1 > 10:
            self.logger.warning("Dramatically different number of resolved nucleotides, using discrepancy -1")
            self.logger.info("Chain %s length %d, chain %s length %d" % (chain1,length1,chain2,length2))
        else:
            # Map units in chain2 to their experimental sequence position ids.
            # These are not experimental sequence positions; different sequences
            # have different ids for the same position
            positions2 = set()           # all positions in chain2
            positions2_to_unit2 = {}     # map those positions back to units
            for unit,position in unit_to_position[chain2].items():
                positions2.add(position)
                positions2_to_unit2[position] = unit

            # Loop over units in chain1, map to positions, and intersect with
            # the positions that go with units in chain2
            # sort by position so it's easier to read and understand
            for unit1,position1 in sorted(unit_to_position[chain1].items(), key=lambda x: x[1]):
                positions1 = set(position_to_position[position1])
                intersection = positions1 & positions2  # intersect positions from 1 and from 2
                positions2 = positions2 - intersection
                if len(intersection) > 1:
                    # could happen because one nucleotide has both A and B alt ids
                    # But I don't see evidence of it happening, which is strange
                    self.logger.info("Trouble: Found multiple matches:")
                    for c in intersection:
                        self.logger.info("Matched %s and %s" % (unit1,positions2_to_unit2[c]))

                elif len(intersection) == 1:
                    for c in intersection:
                        unit2 = positions2_to_unit2[c]
                        matching_pairs.append((unit1,unit2))
                else:
                    self.logger.info("No match for %s position %s to %s" % (unit1,position1,chain2))

            self.logger.info("get_unit_correspondences_intersect: query found %d matching pairs" % len(matching_pairs))

        return matching_pairs


    # def get_unit_correspondences(self, corr_id, info1, info2):
    #     """
    #     Old method, no longer used because it is too slow

    #     Query the database to find all pairs of unit ids with a given
    #     correspondence id, matching two given chains and matching PDB ids,
    #     because that is the fastest query.

    #     Parameters
    #     ----------
    #     corr_id : int
    #         The correspondence id.
    #     info1 : dict
    #         The result of `info` for the first chain to compare.
    #     info2 : dict
    #         The result of `info` for the second chain to compare.

    #     Returns
    #     -------
    #     matching_pairs : list of pairs of unit ids

    #     """

    #     with self.session() as session:
    #         corr_pos = mod.CorrespondencePositions     # C
    #         exp_map1 = aliased(mod.ExpSeqUnitMapping)           # M1
    #         exp_map2 = aliased(mod.ExpSeqUnitMapping)           # M2

    #         mycolumns = [corr_pos.correspondence_id, exp_map1.unit_id.label('unit1'), exp_map2.unit_id.label('unit2')]

    #         # if mycolumns does not reference the first table to be used in the join, use select_from
    #         # select_from(corr_pos).\
    #         # The query without joining units1, units2 takes about 41 seconds on rnatest, for T.th. LSU
    #         # Joining on units1, units2 to get the pdb id takes about 58 seconds
    #         # Not filtering for chain takes many, many minutes, so don't do that
    #         # filtering exp_map1.unit_id.contains(info1['pdb']) takes 46 seconds
    #         # filtering also on info2['pdb'] takes 27 to 39 seconds, so that seems to be the fastest
    #         # filtering with like(info1['pdb']+"%") on both 1 and 2 takes 79 seconds

    #         query = session.query(*mycolumns).\
    #             join(exp_map1, exp_map1.exp_seq_position_id == corr_pos.exp_seq_position_id_1).\
    #             join(exp_map2, exp_map2.exp_seq_position_id == corr_pos.exp_seq_position_id_2).\
    #             filter(corr_pos.correspondence_id == corr_id).\
    #             filter(exp_map1.chain == info1['chain_name']).\
    #             filter(exp_map2.chain == info2['chain_name']).\
    #             filter(exp_map1.unit_id.contains(info1['pdb'])).\
    #             filter(exp_map2.unit_id.contains(info2['pdb']))

    #         """
    #             filter(exp_map1.unit_id.like(info1['pdb']+"%")).\
    #             filter(exp_map2.unit_id.like(info2['pdb']+"%"))

    #             limit(10)
    #             join(units1, units1.unit_id == exp_map1.unit_id).\
    #             join(units2, units2.unit_id == exp_map2.unit_id).\
    #             filter(units1.pdb_id == info1['pdb']).\
    #             filter(units2.pdb_id == info2['pdb'])
    #             filter(units1.unit_id <> units2.unit_id)
    #             filter(units1.chain == info1['chain_name']).\
    #             filter(units2.chain == info2['chain_name'])
    #         """

    #         # the query does not seem to actually run until you ask for data from it;
    #         # setting up is instantaneous, evaluating or counting takes time
    #         self.logger.info("get_unit_correspondences: query set up for %s" % corr_id)

    #         matching_pairs = []
    #         for r in query:
    #             if r.unit1 and r.unit2 and "|" in r.unit1 and "|" in r.unit2:
    #                 matching_pairs.append((r.unit1,r.unit2))

    #         self.logger.info("get_unit_correspondences: query found %d matching pairs" % len(matching_pairs))

    #         return matching_pairs

    def get_unit_correspondences_id_chain(self, corr_id, chain1_name, chain2_name):
        """
        Query the database to find all pairs of unit ids with a given
        correspondence id and matching two given chains, hoping to reduce
        the overall query time across many discrepancy calculations
        This gives more PDB ids than you want at any one time.

        Parameters
        ----------
        corr_id : int
            The correspondence id.
        info1 : dict
            The result of `info` for the first chain to compare.
        info2 : dict
            The result of `info` for the second chain to compare.

        Returns
        -------
        matching_pairs : list of pairs of unit ids

        """

        with self.session() as session:
            corr_pos = mod.CorrespondencePositions     # C
            exp_map1 = aliased(mod.ExpSeqUnitMapping)  # M1
            exp_map2 = aliased(mod.ExpSeqUnitMapping)  # M2

            mycolumns = [corr_pos.correspondence_id, exp_map1.unit_id.label('unit1'), exp_map2.unit_id.label('unit2')]

            # if mycolumns does not reference the first table to be used in the join, use select_from
            # select_from(corr_pos).\
            # The query without joining units1, units2 takes about 41 seconds on rnatest, for T.th. LSU
            # Joining on units1, units2 to get the pdb id takes about 58 seconds
            # Not filtering for chain takes many, many minutes, so don't do that

            query = session.query(*mycolumns).\
                join(exp_map1, exp_map1.exp_seq_position_id == corr_pos.exp_seq_position_id_1).\
                join(exp_map2, exp_map2.exp_seq_position_id == corr_pos.exp_seq_position_id_2).\
                filter(corr_pos.correspondence_id == corr_id).\
                filter(exp_map1.chain == chain1_name).\
                filter(exp_map2.chain == chain2_name)

            # the query does not seem to actually run until you ask for data from it;
            # setting up is instantaneous, evaluating or counting takes time

            matching_pairs = []
            for r in query:
                if r.unit1 and r.unit2 and "|" in r.unit1 and "|" in r.unit2:
                    matching_pairs.append((r.unit1,r.unit2))

            self.logger.info("get_unit_correspondences_id_chain: query found %d matching pairs" % len(matching_pairs))

            return matching_pairs

    def filter_unit_correspondences(self, matching_pairs, info1, info2):
        """
        Filter matching pairs to check that they have the desired
        model and symmetry operator and alt_id.
        Not sure how often the symmetry and alt_id are used.

        This type of checking is superseded by changes on 6/16/2024
        to how the unit_to_position mapping is made.
        This logic does not look very good.
        """

        OK_pairs = []
        for (unit1,unit2) in matching_pairs:
            fields1 = unit1.split("|")
            fields2 = unit2.split("|")

            if fields1[1] == str(info1['model']) and fields2[1] == str(info2['model']):
                if len(fields1) < 7 or fields1[6] == info1['alt_id'] or (fields1[6] == "" and info1['alt_id'] is None):
                    if len(fields2) < 7 or fields2[6] == info2['alt_id'] or (fields2[6] == "" and info2['alt_id'] is None):
                        if len(fields1) < 9 or fields1[8] == info1['sym_op']:
                            if len(fields2) < 9 or fields2[8] == info2['sym_op']:
                                OK_pairs.append((unit1,unit2))

        return OK_pairs


    def load_centers_rotations_pickle(self,info,allunitdictionary):

        for chain_string in info['ife_id'].split('+'):
            # unit directory relative to hub-core
            unit_directory = 'data/units'
            picklefile = os.path.join(unit_directory,chain_string.replace('|','-') + '_NA.pickle')

            if not chain_string in allunitdictionary:

                try:
                    unit_ids, chainIndices, centers, rotations = pickle.load(open(picklefile,"rb"))

                    allunitdictionary[chain_string] = True    # note that this chain was read

                    # add center, rotation pairs to the dictionary according to unit id
                    for i in range(0,len(unit_ids)):
                        if len(centers[i]) == 3 and len(rotations[i]) == 3:
                            allunitdictionary[unit_ids[i]] = (centers[i],rotations[i])
                        else:
                            self.logger.info("Trouble with center/rotation for %s" % unit_ids[i])
                            self.logger.info(str(centers[i]))
                            self.logger.info(str(rotations[i]))

                    self.logger.info('load_centers_rotations_pickle: Loaded %s' % chain_string)

                except:
                    self.logger.info("Could not read pickle file %s " % picklefile)

        return allunitdictionary

    def compare_list_of_pairs(self,old_list,new_list):
        """
        Compare lists of pairs of unit ids from different methods,
        which may put the lists in different orders.
        """

        new_missing = set(old_list) - set(new_list)

        old_missing = set(new_list) - set(old_list)

        if len(old_missing) > 0:
            self.logger.warning('Problem: New list found %s that is not in old list' % str(old_missing))

            for om in old_missing:
                for pp in old_list:
                    if om[0] == pp[0] or om[1] == pp[1]:
                        self.logger.info("This match is in the old list: %s with %s" % pp)

        if len(new_missing) > 0:
            self.logger.warning('Problem: Old list found %s that is not in new list' % str(new_missing))

            for nm in new_missing:
                for pp in new_list:
                    if nm[0] == pp[0] or nm[1] == pp[1]:
                        self.logger.info("This match is in the new list: %s with %s" % pp)

        if len(new_missing) == 0 and len(old_missing) == 0:
            self.logger.info('The old and new methods give identical lists!')


    def discrepancy(self, centers1, centers2, rot1, rot2):
        """
        Compute the discrepancy using centers and rotation matrices.

        Parameters
        ----------
        centers1 : list
            List of numpy.array of the centers for first chain.
        centers2 : list
            List of numpy.array of the centers for second chain.
        rot1 : list
            List of numpy.array of the rotation matrices for first chain.
        rot2 : list
            List of numpy.array of the rotation matrices for second chain.

        Returns
        -------
        data : float
            The computed discrepancy.
        """

        disc = matrix_discrepancy(centers1, rot1, centers2, rot2)

        if np.isnan(disc):
            raise core.InvalidState("NaN for discrepancy")

        return disc


    def get_chain_info(self, chain_id):
        """Load the required information about a chain. Since we want to use
        the results of this loader for the NR stages we use the same data as
        was in the IFE's the given chain is a part of.

        Parameters
        ----------
        chain_id : int
            The chain id to look up.

        Returns
        -------
        ife_info : dict
            A dict with a 'chain_name', 'chain_id', `pdb`, `model`, `ife_id`,
            `sym_op`, 'alt_id', and `name` keys.
        """

        with self.session() as session:
            query = session.query(mod.ChainInfo.chain_name,
                                  mod.ChainInfo.chain_id,
                                  mod.ChainInfo.chain_length,
                                  mod.IfeInfo.pdb_id.label('pdb'),
                                  mod.IfeInfo.model,
                                  mod.IfeInfo.ife_id,
                                  ).\
                join(mod.IfeChains,
                     mod.IfeChains.chain_id == mod.ChainInfo.chain_id).\
                join(mod.IfeInfo,
                     mod.IfeInfo.ife_id == mod.IfeChains.ife_id).\
                filter(mod.IfeInfo.new_style == True).\
                filter(mod.ChainInfo.chain_id == chain_id)

            # the following lines were not indented this much, changed 2022-05-18
            result = [row for row in query]

            #if not query.count():
            if len(result) == 0:
                raise core.InvalidState("Could not load chain with id %s" %
                                        chain_id)
            ife = ut.row2dict(query.first())

        with self.session() as session:
            query = session.query(mod.UnitInfo.sym_op,
                                  mod.UnitInfo.alt_id).\
                join(mod.ChainInfo,
                     (mod.ChainInfo.pdb_id == mod.UnitInfo.pdb_id) &
                     (mod.ChainInfo.chain_name == mod.UnitInfo.chain)).\
                join(mod.UnitCenters,
                     mod.UnitCenters.unit_id == mod.UnitInfo.unit_id).\
                join(mod.UnitRotations,
                     mod.UnitRotations.unit_id == mod.UnitInfo.unit_id).\
                filter(mod.ChainInfo.chain_id == chain_id).\
                distinct()

            result = [row for row in query]

            #if not query.count():       # caused QueuePool limit of size 40 error
            if len(result) == 0:
                raise core.InvalidState("Could not get info for chain %s" %
                                        chain_id)

            ife['sym_op'] = pick(['1_555', 'P_1'], 'sym_op', query)
            ife['alt_id'] = pick([None, 'A', 'B'], 'alt_id', query)
            ife['name'] = ife['ife_id'] + '+' + ife['sym_op']
            return ife


    def __correspondence_query__(self, chain1, chain2):
        """Create a query for correspondences between the two chains. This only
        checks in the given direction.

        Parameters
        ----------
        chain1 : int
            The first chain id.
        chain2 : int
            The second chain id.

        Returns
        -------
        corr_id : int
            The correspondence id if there is an alignment between the two
            chains.
        """

        with self.session() as session:
            info = mod.CorrespondenceInfo
            mapping1 = aliased(mod.ExpSeqChainMapping)
            mapping2 = aliased(mod.ExpSeqChainMapping)
            query = session.query(info.correspondence_id).\
                join(mapping1, mapping1.exp_seq_id == info.exp_seq_id_1).\
                join(mapping2, mapping2.exp_seq_id == info.exp_seq_id_2).\
                filter(mapping1.chain_id == chain1).\
                filter(mapping2.chain_id == chain2).\
                filter(info.good_alignment >= 1)

            result = query.first()
            if result:
                return result.correspondence_id
            return None


    def corr_id(self, chain_id1, chain_id2):
        """Given two chain ids, load the correspondence id between them. This
        will raise an exception if these chains have not been aligned, or if
        there is more than one alignment between these chains. The chain ids
        should be the ids used in the database.

        Parameters
        ----------
        chain_id1 : int
            First chain id.
        chain_id2 : int
            Second chain id.

        Raises
        ------
        core.Skip
            If there is no alignment between the two chains

        Returns
        -------
        corr_id : int
            The correspondence id for the alignment between the two chains.
        """

        self.logger.debug("corr_id: chain_id1: %s" % chain_id1)
        self.logger.debug("corr_id: chain_id2: %s" % chain_id2)

        corr_id = self.__correspondence_query__(chain_id1, chain_id2)
        if corr_id is None:
            corr_id = self.__correspondence_query__(chain_id2, chain_id1)

        return corr_id


    def __check_matrices__(self, table, info):
        """Check the that there are entries for the given info dict in the
        given table.

        Parameters
        ----------
        table : Table
            The table to look up in
        info : dict
            The info dict to use.

        Returns
        -------
        has_entries : bool
            True if there are entries in the table.
        """

        with self.session() as session:
            query = session.query(table).\
                join(mod.UnitInfo,
                     mod.UnitInfo.unit_id == table.unit_id).\
                filter(mod.UnitInfo.pdb_id == info['pdb']).\
                filter(mod.UnitInfo.chain == info['chain_name']).\
                filter(mod.UnitInfo.model == info['model']).\
                filter(mod.UnitInfo.sym_op == info['sym_op']).\
                limit(1)

            return bool(query.count())


    def has_matrices(self, info):
        """Check if the given info has all the required matrices. This will
        check in the database for entries in the unit_centers and
        unit_rotations for the given chain.

        Parameters
        ----------
        info : dict
            The info dict to check.

        Returns
        -------
        valid : bool
            True if the given chain has all required matrices.
        """

        return self.__check_matrices__(mod.UnitCenters, info) and \
            self.__check_matrices__(mod.UnitRotations, info)

    def gather_matching_centers_rotations(self,unit_pairs,allunitdictionary):

        """
        loop over pairs of unit ids, pull out corresponding center and
        rotation matrices for the pairs of unit ids.
        """

        c1 = []
        c2 = []
        r1 = []
        r2 = []
        seen1 = set()
        seen2 = set()

        for (unit1,unit2) in unit_pairs:

            if unit1 in seen1:
                #raise core.InvalidState("gather_matching_centers_rotations: Got duplicate unit1 %s" % unit1)
                self.logger.info("gather_matching_centers_rotations: Got duplicate unit1 %s" % unit1)

            elif unit2 in seen2:
                #raise core.InvalidState("gather_matching_centers_rotations: Got duplicate unit2 %s" % unit2)
                self.logger.info("gather_matching_centers_rotations: Got duplicate unit2 %s" % unit2)

            else:
                # only add when neither has been seen already
                seen1.add(unit1)
                seen2.add(unit2)

                if allunitdictionary.get(unit1) is not None and allunitdictionary.get(unit2) is not None:
                    c1.append(allunitdictionary[unit1][0])
                    c2.append(allunitdictionary[unit2][0])
                    r1.append(allunitdictionary[unit1][1])
                    r2.append(allunitdictionary[unit2][1])

        return np.array(c1), np.array(c2), np.array(r1), np.array(r2)


    def calculate_discrepancy(self, info1, info2, corr_id, c1, c2, r1, r2):
        """
        Compute the discrepancy between two given chains.
        This will produce a list of two dictionaries containing data for
        the database, in both orders.

        Parameters
        ----------
        info1 : dict
            The first chain to use.
        info2 : dict
            The second chain to use
        corr_id : dict
            The correspondence id.
        c1, c2: list of 3-dimensional vectors
            Centers of matching nucleotides
        r1, r2: list of 3x3 rotation matrices
            Rotation matrices of matching nucleotides

        Returns
        -------
        entries : list
            A list of dictonaries with the discrepancies betweeen these chains.
            The dictonaries will have keys `chain_id_1`, `chain_id_2`,
            `model_1`, `model_2`, `correspondence_id`, `discrepancy` and
            `num_nucleotides`.
        """

        if len(c1) == 0:
            self.logger.warning("No matched nucleotides, using discrepancy = -1")
            disc = -1
        elif len(c1) < 3:
            self.logger.warning("Too few matched nucleotides to compute discrepancy for %s %s, using disc = -1 instead" %
                                  (info1['name'], info2['name']))
            disc = -1
        else:

            try:
                disc = matrix_discrepancy(c1, r1, c2, r2)

            except Exception as err:
                self.logger.warning("Could not compute discrepancy for %s %s, using discrepancy = -1 instead" %
                                  (info1['name'], info2['name']))
                disc = -1

            if np.isnan(disc):
                self.logger.warning("Could not compute discrepancy for %s %s, using discrepancy = -1 instead" %
                                  (info1['name'], info2['name']))
                disc = -1


        entry1 = {
            'chain_id_1': info1['chain_id'],
            'chain_id_2': info2['chain_id'],
            'model_1': info1['model'],
            'model_2': info2['model'],
            'correspondence_id': corr_id,
            'discrepancy': float(disc),
            'num_nucleotides': len(c1)
        }

        entry2 = {
            'chain_id_1': info2['chain_id'],
            'chain_id_2': info1['chain_id'],
            'model_1': info2['model'],
            'model_2': info1['model'],
            'correspondence_id': corr_id,
            'discrepancy': float(disc),
            'num_nucleotides': len(c1)
        }

        return [entry1, entry2]

    def data(self, chain_ids, **kwargs):
        """
        New in December 2020.
        Compute all chain to chain discrepancies in the given group.
        Loop over all pairs of chain ids in this group.

        This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.

        Parameters
        ----------
        chain_ids : list of integer chain ids
            The list of chain ids to compute data for

        Returns
        -------
        data : list
            The discrepancies of each pair of chains.
        """

        L = len(chain_ids)
        self.logger.info("Computing discrepancies for a group of %d chains" % L)
        self.logger.info("%d discrepancies needed in this group" % (L*(L-1)/2))

        # from the list of chain_ids, find all pairs that already have a discrepancy computed
        already_computed = []
        already_computed_discrepancy = []
        sim = mod.ChainChainSimilarity
        with self.session() as session:
            query = session.query(sim).\
                filter(or_(and_(sim.chain_id_1.in_(chain_ids),sim.chain_id_2.in_(chain_ids))))

            for r in query:
                if r.chain_id_1 != r.chain_id_2:
                    already_computed.append((r.chain_id_1,r.chain_id_2))
                    already_computed_discrepancy.append((r.chain_id_1,r.chain_id_2,r.discrepancy))
        self.logger.info("Found %d discrepancy values already calculated" % (len(already_computed)/2))

        # debugging
        # figure out why more discrepancies are computed than are needed in some cases
        if len(already_computed)/2 > L*(L-1)/2:
            self.logger.info("Many groups on production have duplicate entries of discrepancies.  Not sure why.  Sometimes the discrepancy differs as well.")

            # uncomment the code below to see many, many examples
            if 0 > 1:
                already_computed_discrepancy = sorted(already_computed_discrepancy)
                for i in range(0,len(already_computed_discrepancy)-1):
                    if already_computed_discrepancy[i][1] == already_computed_discrepancy[i+1][1] and already_computed_discrepancy[i][0] == already_computed_discrepancy[i+1][0]:
                        self.logger.info("Duplicate entry of chains %d and %d with discrepancies %10.7f and %10.7f" % (already_computed_discrepancy[i][0], already_computed_discrepancy[i][1], already_computed_discrepancy[i][2], already_computed_discrepancy[i+1][2]))

        # loop over pairs of chain ids in this group, skipping those that already have discrepancy calculated, look up corr_id
        check_pairs = sorted(list(set(it.combinations(chain_ids, 2)) - set(already_computed)))
        if len(check_pairs) > 0:
            self.logger.info("Looking up correspondence ids for %d pairs of chains that need to be computed" % len(check_pairs))

        required_pairs = []
        log_count = 0
        chain_info = {}

        for (chain1_id,chain2_id) in check_pairs:
            corr_id = self.corr_id(chain1_id, chain2_id)
            if corr_id is None:
                if log_count < 20:
                    if not chain1_id in chain_info:
                        chain_info[chain1_id] = self.get_chain_info(chain1_id)
                    if not chain2_id in chain_info:
                        chain_info[chain2_id] = self.get_chain_info(chain2_id)
                    info1 = chain_info[chain1_id]
                    info2 = chain_info[chain2_id]

                    self.logger.info("No correspondence id between chains %s %s and %s %s" % (info1['ife_id'],chain1_id,info2['ife_id'],chain2_id))
                    log_count = log_count + 1
                else:
                    self.logger.info("No correspondence id between chains %s and %s" % (chain1_id,chain2_id))
                # Note: cannot store a discrepancy with a null correspondence id
            else:
                required_pairs.append((corr_id,chain1_id,chain2_id))
        self.logger.info("Found %d discrepancy values needing to be calculated" % len(required_pairs))



        # This block is for debugging, to check to see if we get the same discrepancies as before
        Recompute = True   # run the debugging
        Recompute = False

        # check to see that we get the same discrepancy as before for some cases
        # this is for debugging; generally the program will be run with Recompute = False
        if Recompute and len(already_computed_discrepancy) > 0:

            # retrieve chain information for all chains once and store the data
            # 50 seconds for T.th. SSU on rnatest in December 2020

            self.logger.info("Retrieving chain information once for each chain")
            chain_info = {}
            for chain_id in chain_ids:
                chain_info[chain_id] = self.get_chain_info(chain_id)

            # takes about 44 seconds on rnatest in December 2020
            # takes about 2 seconds on rnaprod2 in September 2024
            self.logger.info("Loading unit to position correspondences")
            unit_to_position = pickle.load(open("unit_to_position.pickle","rb"))

            # takes about 45 seconds on rnatest in December 2020
            # takes about 7 seconds on rnaprod2 in September 2024
            self.logger.info("Loading position to position correspondences")
            position_to_position = pickle.load(open("position_to_position.pickle","rb"))

            # Loop over needed pairs of chains, query for unit correspondences, and calculate discrepancies
            # The slowest part of the process is the query for unit correspondences.
            # One hope was that querying by corr_id and the two chains and then narrowing down to the desired
            # PDBs would be faster than doing them individually, but that is sometimes much slower.

            allunitdictionary = defaultdict()            # store up centers and rotations

            current = 1
            chain1_seen = set()

            LL = min(len(already_computed_discrepancy),50)

            for (chain1_id,chain2_id,discrepancy) in already_computed_discrepancy[0:LL]:

                # make sure not to repeat any pairings in different order
                if chain1_id < chain2_id:

                    self.logger.info("Re-computing discrepancy %d for this group" % (current))
                    current += 1

                    chain1_seen.add(chain1_id)
                    if len(chain1_seen) > 20:
                        chain1_seen = set()
                        allunitdictionary = defaultdict()    # avoid accumulating data forever; reset sometimes

                    info1 = chain_info[chain1_id]
                    info2 = chain_info[chain2_id]

                    # new method
                    self.logger.info("Intersect for matching units for chain %s, chain %s" % (info1['ife_id'],info2['ife_id']))
                    unit_pairs = self.get_unit_correspondences_intersect(info1,info2,unit_to_position,position_to_position)

                    # filter out units with wrong symmetry or alt id
                    # unit_pairs = self.filter_unit_correspondences(new_unit_pairs,info1,info2)

                    # show some matched units to build confidence
                    if len(unit_pairs) > 0:
                        for i in range(0,min(5,len(unit_pairs))):
                            self.logger.info("Matched unit ids %s and %s" % unit_pairs[i])

                        # load center and rotation data for the current ifes, if not already loaded
                        allunitdictionary = self.load_centers_rotations_pickle(info1,allunitdictionary)
                        allunitdictionary = self.load_centers_rotations_pickle(info2,allunitdictionary)


                    # gather matching centers and rotations for these chains
                    [c1, c2, r1, r2] = self.gather_matching_centers_rotations(unit_pairs,allunitdictionary)
                    self.logger.info("Got matching centers and rotations")

                    # compute the discrepancy between these IFEs
                    # if wrong numbers of matched nucleotides, discrepancy will be -1
                    corr_id = 0
                    discrepancies = self.calculate_discrepancy(info1, info2, corr_id, c1, c2, r1, r2)

                    self.logger.info("Previous discrepancy %0.7f" % (discrepancy))
                    self.logger.info("New      discrepancy %0.7f" % (discrepancies[0]["discrepancy"]))

                    if abs(10000000*(discrepancy-discrepancies[0]["discrepancy"])) > 1000:
                        self.logger.info("Problem: old and new discrepancies really don't agree")

                    if abs(10000000*(discrepancy-discrepancies[0]["discrepancy"])) > 10:
                        self.logger.info("Problem: old and new discrepancies don't agree")




        if not Recompute and len(required_pairs) > 0:

            # takes about 44 seconds on rnatest in December 2020
            self.logger.info("Loading unit to position correspondences")
            unit_to_position = pickle.load(open("unit_to_position.pickle","rb"))

            # takes about 45 seconds on rnatest in December 2020
            self.logger.info("Loading position to position correspondences")
            position_to_position = pickle.load(open("position_to_position.pickle","rb"))

            # Loop over needed pairs of chains, query for unit correspondences, and calculate discrepancies
            # The slowest part of the process is the query for unit correspondences.
            # One hope was that querying by corr_id and the two chains and then narrowing down to the desired
            # PDBs would be faster than doing them individually, but that is sometimes much slower.

            allunitdictionary = defaultdict()            # store up centers and rotations

            current = 0
            chain1_seen = set()                          # count how many times a specific chain1 is seen
            for (corr_id,chain1_id,chain2_id) in required_pairs:

                current += 1
                self.logger.info("Computing discrepancy %d of %d for this group" % (current,len(required_pairs)))

                chain1_seen.add(chain1_id)
                if len(chain1_seen) > 20:
                    chain1_seen = set()
                    allunitdictionary = defaultdict()    # avoid accumulating data forever; reset sometimes

                # Store chain info so that each one is only looked up once.
                # It would seem to be better to just retrieve all of them at once,
                # but that crashed the pipeline once with "QueuePool limit overflow"
                # so it seems better to do them one by one.

                if not chain1_id in chain_info:
                    chain_info[chain1_id] = self.get_chain_info(chain1_id)
                if not chain2_id in chain_info:
                    chain_info[chain2_id] = self.get_chain_info(chain2_id)

                info1 = chain_info[chain1_id]
                info2 = chain_info[chain2_id]

                # Future work:
                # Recognize multiple chains for IFEs made of more than one chain; currently only 1st chain is used

                # new method
                self.logger.info("Intersect for matching units for chain %s, chain %s" % (info1['ife_id'],info2['ife_id']))
                unit_pairs = self.get_unit_correspondences_intersect(info1,info2,unit_to_position,position_to_position)

                # filter out units with wrong symmetry or alt id
                unit_pairs = self.filter_unit_correspondences(unit_pairs,info1,info2)

                # Next: if matched resolved nucleotides is much less than one of the chain lengths, discrepancy = NaN

                """
                # old method for getting correspondences
                self.logger.info("Query for matching units for chain %s, chain %s" % (info1['ife_id'],info2['ife_id']))
                old_unit_pairs = self.get_unit_correspondences(corr_id,info1,info2)

                # filter out units with wrong symmetry or alt id
                old_unit_pairs = self.filter_unit_correspondences(old_unit_pairs,info1,info2)
                self.logger.info("Found %d unit id pairs for chain %s chain %s" % (len(unit_pairs),info1['ife_id'],info2['ife_id']))

                self.compare_list_of_pairs(old_unit_pairs,unit_pairs)
                """

                # show some matched units to build confidence
                if len(unit_pairs) > 0:
                    for i in range(0,min(5,len(unit_pairs))):
                        self.logger.info("Matched unit ids %s and %s" % unit_pairs[i])

                    # load center and rotation data for the current ifes, if not already loaded
                    allunitdictionary = self.load_centers_rotations_pickle(info1,allunitdictionary)
                    allunitdictionary = self.load_centers_rotations_pickle(info2,allunitdictionary)

                # gather matching centers and rotations for these chains
                [c1, c2, r1, r2] = self.gather_matching_centers_rotations(unit_pairs,allunitdictionary)
                self.logger.info("Gathered %d matching centers and rotations for %s and %s, %d of %d in this group" % (len(c1),info1['ife_id'],info2['ife_id'],current,len(required_pairs)))

                # compute the discrepancy between these IFEs
                # if wrong numbers of matched nucleotides, discrepancy will be -1
                discrepancies = self.calculate_discrepancy(info1, info2, corr_id, c1, c2, r1, r2)

                self.logger.info("Discrepancy to load: %s" % discrepancies[0])
                for d in discrepancies:
                    yield mod.ChainChainSimilarity(**d)
