"""
Compute chain to chain discrepancies within equivalence classes.
This will 1) look at good
correspondences, 2) extract all the aligned chains, 3) compute the
geometric discrepancy between them, and 4) place the discrepanices
in the database.

This will only compare the first chain in each IFE.
"""

import functools as ft
import itertools as it
import numpy as np
import operator as op
import pickle
import time

from collections import defaultdict

from sqlalchemy import and_
from sqlalchemy import or_
from sqlalchemy.orm import aliased
from sqlalchemy.sql import union_all

from fr3d.geometry.discrepancy import matrix_discrepancy

from pymotifs import core
import pymotifs.utils as ut
from pymotifs import models as mod
from pymotifs.utils import discrepancy as disc

from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqUnitMappingLoader
from pymotifs.ife.loader import Loader as IfeLoader
from pymotifs.units.centers import Loader as CenterLoader
from pymotifs.units.rotation import Loader as RotationLoader

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


def label_center(table, number):
    """Select the center columns from the given table. This will alias the
    columns to ones that are unique

    Parameters
    ----------
    table : Table
        The table to select columns from.
    number : int
        The number to postfix

    Returns
    -------
    columns : list
        The list of columns to select.
    """

    return [
        getattr(table, 'unit_id').label('unit%s' % number),
        getattr(table, 'x').label('x%s' % number),
        getattr(table, 'y').label('y%s' % number),
        getattr(table, 'z').label('z%s' % number),
    ]



def make_unique_list(input_list):
    unique_list = []
    for foo in input_list:
        if foo not in unique_list:
            unique_list.append(foo)
    return unique_list


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
                        IfeLoader, CenterLoader, RotationLoader])


    def known_unit_entries(self, table):
        """Create a set of (pdb, chain) tuples for all chains that have entries
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
        with self.session() as session:
            info = mod.UnitInfo
            query = session.query(info.pdb_id,
                                  info.chain,
                                  ).\
                join(table, table.unit_id == info.unit_id).\
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
        """This will compute all pairs to compare. This will group all pdbs
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


        GeneratePickleFiles = False   # appropriate to use when debugging the rest of the program
        GeneratePickleFiles = True    # must be used in production, to update the files each week

        if GeneratePickleFiles:
            # query and write to disk unit correspondence data once per run of chain_chain/comparison.py
            # takes about 2 1/3 minutes on production in December 2020
            self.logger.info('Getting the unit_id to experimental sequence position table')
            with self.session() as session:
                EM = mod.ExpSeqUnitMapping
                query = session.query(EM.unit_id,EM.exp_seq_position_id).select_from(EM)

                unit_to_position = {}      # will be a dictionary of dictionaries
                count = 0
                for r in query:
                    if r.unit_id and "|" in r.unit_id:    # not sure why, but some rows have None
                        fields = r.unit_id.split("|")
                        if len(fields) > 3:
                            key = "|".join(fields[0:3])   # pdb id, model, chain is the top level key
                            if not key in unit_to_position:
                                unit_to_position[key] = {}
                            unit_to_position[key][r.unit_id] = r.exp_seq_position_id
                            count += 1
                self.logger.info('Got %d unit to position mappings' % count)
                pickle.dump(unit_to_position, open("unit_to_position.pickle", "wb" ), 2)
                self.logger.info('Wrote pickle file unit_to_position.pickle')

            unit_to_position = {}     # hopefully clear this from memory

        if GeneratePickleFiles:
            # takes about 6.5 minutes on production in December 2020
            self.logger.info('Getting the position to position mappings')

            position_to_position = defaultdict(list)
            count = 0

            i = 0
            j = 1
            newcount = 1

            while j < 9999:
                if newcount == 0:
                    j = 9999            # this will be the last query, check for very large numbers
                else:
                    j = i+1

                with self.session() as session:
                    CP = mod.CorrespondencePositions
                    query = session.query(CP.exp_seq_position_id_1,CP.exp_seq_position_id_2).\
                    select_from(CP).\
                    filter(CP.exp_seq_position_id_1 < 100000*j).\
                    filter(CP.exp_seq_position_id_1 >= 100000*i)

                    newcount = 0
                    for r in query:
                        position_to_position[r.exp_seq_position_id_1].append(r.exp_seq_position_id_2)
                        count += 1
                        newcount += 1

                    self.logger.info('Total of %10d position to position correspondences, id up to %d' % (count,100000*j))

                i = i + 1

            pickle.dump(position_to_position, open("position_to_position.pickle", "wb" ), 2)
            self.logger.info('Wrote pickle file position_to_position.pickle')
            self.logger.info('Maximum exp_seq_position_id_1 is %d' % max(position_to_position.keys()))

            position_to_position = {}   # hopefully clear the memory


        # to speed things up when debugging, return this (first,seconds) list for T.th. LSU
        # self.logger.info("Skipping grouping, just using a small T.th. LSU group")
        # return [(92730, [91729, 91730, 91837, 91838, 67402, 43569, 67346, 43811, 43755, 43866, 61884, 61829, 87224, 44080, 44023, 44136, 91676, 45731, 45790, 53900, 43432, 50768, 43486, 45098, 45042, 57600, 83493, 46559, 57539, 60287, 46795, 46500, 59606, 46736, 56239, 50710, 60548, 46139, 56178, 60495, 59498, 46196, 60296, 90105, 37220, 64813, 84743, 58749, 37783, 84235, 64239, 86194, 46677, 84684, 84257, 46618, 59323, 59215, 37275, 54309, 92562, 46253, 54117, 56300, 45668, 37838, 45606, 46382, 62569, 59107, 59388, 54252, 62553, 55547, 92618, 64921, 46441, 87643, 64705, 56483, 87616, 47757, 56603, 90212, 55655, 59866, 65030, 51907, 92674, 51958, 56422, 59507, 64822, 56361, 56543, 53211, 65039, 90222, 59875, 37108, 55494, 59615, 47586])]

#        self.logger.info("Skipping grouping on rnatest, just using 20 chains from the largest T.th. SSU group")
#        return [[60295, 53866, 50732, 50674, 43488, 43721, 45008, 43833, 46045, 45064, 64821, 43777, 45754, 43434, 43992, 45695, 1766, 44105, 43594, 44048, 1791, 61917, 61862, 83490, 87209, 67344, 876, 91677, 67400, 429, 91839, 91840, 91731, 91732]]

#        self.logger.info("Skipping grouping on rnatest, just using the largest T.th. SSU group")
#        return [[27075, 1002, 68491, 928, 70901, 1720, 724, 24733, 70846, 656, 27602, 27581, 63337, 68342, 58018, 24756, 65833, 70985, 678, 58074, 67255, 27560, 63360, 25324, 55970, 60698, 71095, 70736, 61541, 24802, 60751, 65725, 34161, 27623, 57796, 65779, 27518, 67143, 906, 27665, 57684, 8284, 18221, 25450, 34273, 57659, 85003, 71040, 69477, 63314, 57740, 58494, 24779, 65887, 70791, 71150, 63120, 65251, 61408, 63868, 34329, 84950, 64128, 700, 952, 57852, 20158, 58435, 85056, 77040, 60804, 85109, 64184, 67091, 83763, 87947, 22048, 63602, 23041, 22991, 20893, 4152, 63383, 88031, 83935, 83823, 88003, 4175, 4198, 20805, 66481, 61714, 61660, 53444, 38009, 67199, 61596, 36850, 38345, 83711, 66651, 38400, 23014, 85479, 27686, 63487, 66706, 66423, 23199, 37298, 38457, 37522, 77096, 84016, 66926, 47437, 20870, 37635, 38512, 61464, 57963, 46582, 65306, 66596, 57908, 52958, 83879, 47493, 67311, 5537, 47107, 85425, 53497, 84102, 46641, 11002, 51979, 66816, 47163, 84838, 59023, 21880, 37353, 66761, 83991, 20826, 88100, 52032, 25366, 88077, 84894, 27644, 57340, 46405, 23222, 86054, 38233, 60860, 23176, 84148, 6135, 34217, 59734, 38064, 34385, 82270, 58663, 37577, 38288, 36905, 37691, 66541, 66871, 48157, 52738, 52878, 25429, 57395, 46346, 57450, 52306, 25408, 4129, 25471, 38121, 58531, 22069, 8307, 85943, 18198, 59756, 34495, 54759, 63064, 56028, 48132, 38176, 37465, 50848, 46759, 88123, 69422, 46464, 50906, 46700, 50992, 58905, 63547, 90006, 52903, 71773, 37074, 1743, 54814, 57261, 88054, 86108, 57207, 54274, 57285, 34440, 68210, 52281, 92473, 46103, 45225, 87725, 50790, 87742, 8261, 37410, 36962, 52853, 67984, 56144, 56086, 45282, 53098, 54217, 89997, 37017, 84125, 52661, 34551, 52683, 45571, 45633, 854, 58522, 92419, 56449, 71845, 68201, 59080, 63431, 59874, 71821, 22006, 52640, 56205, 53178, 58721, 27193, 20784, 53043, 54621, 34051, 20847, 53607, 27539, 53550, 84039, 92365, 48216, 48270, 64347, 55623, 56570, 50475, 82261, 56388, 46162, 46219, 44324, 45396, 37129, 58964, 55463, 50533, 84081, 84060, 85934, 34105, 56510, 47723, 45869, 45339, 812, 68093, 92695, 64356, 54868, 44268, 62549, 62565, 55569, 56266, 47550, 45928, 37186, 64024, 83645, 65038, 20918, 56327, 59223, 21922, 27172, 83589, 55516, 56685, 25345, 37749, 59387, 56630, 71969, 46523, 25303, 71797, 92311, 54331, 54925, 71920, 65029, 64076, 64929, 92639, 71870, 59396, 59506, 64920, 60494, 54141, 54390, 90211, 59115, 59214, 64704, 53285, 58757, 27250, 51571, 92527, 37804, 21985, 50591, 68102, 71945, 50649, 71895, 59106, 22027, 64238, 59865, 51622, 59327, 92583, 67993, 58748, 59322, 90221, 51673, 45453, 59497, 64713, 90113, 21901, 45512, 51724, 53153, 53338, 51826, 5610, 90104, 21943, 52138, 833, 57505, 51775, 52085, 60547, 53391, 64247, 37241, 86162, 86215, 47667, 53750, 25387, 59605, 53808, 64812, 45987, 60286, 87611, 87638, 57566, 84649, 51877, 84708, 59614, 84230, 84252, 17109, 51928, 21964, 53232, 60295, 53866, 50732, 50674, 43488, 43721, 45008, 43833, 46045, 45064, 64821, 43777, 45754, 43434, 43992, 45695, 1766, 44105, 43594, 44048, 1791, 61917, 61862, 83490, 87209, 67344, 876, 91677, 67400, 429, 91839, 91840, 91731, 91732]]

#        self.logger.info("Skipping grouping on rnatest, just using a group of 350 chains")
#        return [[57418L, 57308L, 71006L, 38366L, 71116L, 38030L, 80391L, 76219L, 37825L, 73358L, 36300L, 50042L, 45592L, 57473L, 57363L, 62808L, 50813L, 47130L, 47460L, 52283L, 53180L, 56107L, 47745L, 47186L, 47516L, 52880L, 92548L, 63568L, 92440L, 38142L, 92660L, 80271L, 38197L, 92604L, 92386L, 38421L, 71061L, 92494L, 71171L, 92716L, 73483L, 74619L, 74735L, 74793L, 72414L, 45654L, 61431L, 64207L, 70520L, 70637L, 70579L, 70521L, 52685L, 47131L, 47461L, 53045L, 54162L, 48272L, 52740L, 52855L, 53628L, 47517L, 47187L, 88513L, 53100L, 62482L, 62487L, 88535L, 92332L, 38254L, 69443L, 37770L, 37543L, 76159L, 80511L, 80332L, 84529L, 70812L, 37374L, 63141L, 37712L, 37598L, 70922L, 63623L, 90791L, 90682L, 90792L, 90500L, 66783L, 66563L, 66673L, 60827L, 59337L, 83902L, 77119L, 60883L, 88026L, 61618L, 12022L, 57819L, 57930L, 57707L, 34296L, 58041L, 58097L, 34352L, 57875L, 57763L, 90126L, 74677L, 72356L, 72298L, 44106L, 45246L, 64151L, 68224L, 61487L, 73592L, 34408L, 34184L, 34518L, 34463L, 58517L, 70696L, 52905L, 56049L, 47689L, 52308L, 50871L, 26577L, 34974L, 34902L, 63452L, 37319L, 70867L, 37656L, 84528L, 80331L, 80272L, 80512L, 80572L, 80452L, 38309L, 69498L, 38085L, 63508L, 90442L, 84647L, 83611L, 76091L, 36378L, 66444L, 87970L, 61563L, 67166L, 61681L, 66893L, 77063L, 68006L, 90235L, 64726L, 58770L, 59627L, 59409L, 68223L, 64942L, 65051L, 66502L, 61735L, 67334L, 66948L, 66618L, 66838L, 67222L, 66728L, 12023L, 49296L, 49576L, 48821L, 59045L, 73519L, 72529L, 82922L, 70578L, 70697L, 52960L, 54946L, 34986L, 36983L, 70757L, 37431L, 37095L, 34072L, 63085L, 37207L, 80571L, 80451L, 37262L, 37150L, 37486L, 37038L, 34126L, 84646L, 76090L, 90382L, 59128L, 58458L, 48271L, 83958L, 83846L, 65273L, 67278L, 64834L, 58544L, 64260L, 68115L, 90019L, 58545L, 59410L, 64943L, 65328L, 84014L, 49352L, 48647L, 48705L, 57985L, 59102L, 71993L, 47630L, 44289L, 12197L, 44345L, 45655L, 59237L, 60308L, 72243L, 86652L, 86534L, 86475L, 86357L, 48217L, 46008L, 54889L, 50928L, 50697L, 62469L, 51014L, 76220L, 62561L, 83667L, 34240L, 34574L, 60299L, 90127L, 60307L, 68007L, 64261L, 68116L, 90020L, 49063L, 49520L, 49632L, 49464L, 49408L, 48763L, 48912L, 5665L, 5666L, 71919L, 71968L, 45303L, 59236L, 86416L, 35563L, 86593L, 62030L, 70638L, 53155L, 48273L, 76160L, 80392L, 87080L, 87194L, 64727L, 90236L, 59338L, 67036L, 49856L, 49800L, 49912L, 71869L, 43993L, 69408L, 53571L, 46066L, 19153L, 64835L, 58771L, 59628L, 65052L, 49688L, 92605L, 87652L, 87676L, 74868L, 74894L, 87196L, 63453L, 63086L, 63142L, 63509L, 92717L, 59129L, 49744L, 92549L, 62999L, 62909L, 48218L, 50755L, 92661L, 90558L, 48623L, 45593L, 54163L, 84271L, 84289L, 53887L, 48956L, 76101L, 86594L, 86358L, 86417L, 86535L, 86476L, 52856L, 52686L, 52309L, 52881L, 87626L, 87653L, 52284L, 52741L, 84245L, 84267L, 86653L, 62031L, 87135L, 17719L, 46009L, 46067L]]
#        return [[52741L, 84245L, 84267L, 86653L, 62031L, 87135L, 17719L, 46009L, 46067L]]

        # large set of chains on production
#        return [[5665L, 5666L, 12022L, 12023L, 12197L, 17719L, 19153L, 26577L, 34176L, 34230L, 34288L, 34344L, 34400L, 34456L, 34512L, 34567L, 34622L, 34678L, 35030L, 35102L, 35114L, 35691L, 36434L, 36512L, 36837L, 36938L, 37028L, 37523L, 37578L, 37635L, 37690L, 37747L, 37802L, 37859L, 37914L, 37971L, 38026L, 38083L, 38138L, 38196L, 38252L, 38310L, 38365L, 38578L, 38633L, 38690L, 38745L, 38802L, 38857L, 38914L, 38969L, 44094L, 44207L, 44390L, 44446L, 45347L, 45404L, 45693L, 45694L, 45755L, 45756L, 46109L, 46110L, 46167L, 46168L, 47231L, 47232L, 47287L, 47288L, 47561L, 47562L, 47617L, 47618L, 47731L, 47790L, 47846L, 48318L, 48319L, 48372L, 48373L, 48374L, 48724L, 48748L, 48806L, 48864L, 48922L, 49013L, 49057L, 49164L, 49397L, 49453L, 49509L, 49565L, 49621L, 49677L, 49733L, 49789L, 49845L, 49901L, 49957L, 50013L, 50143L, 50798L, 50856L, 50914L, 50972L, 51029L, 51115L, 52384L, 52385L, 52409L, 52410L, 52786L, 52787L, 52841L, 52842L, 52956L, 52957L, 52981L, 52982L, 53006L, 53061L, 53146L, 53201L, 53256L, 53281L, 53672L, 53729L, 53988L, 54263L, 54264L, 54990L, 55047L, 56150L, 56208L, 57409L, 57464L, 57519L, 57574L, 57808L, 57864L, 57920L, 57976L, 58031L, 58086L, 58142L, 58198L, 58252L, 58253L, 58308L, 58309L, 58809L, 58868L, 58895L, 58896L, 59121L, 59122L, 59396L, 59453L, 59479L, 59480L, 59587L, 59588L, 59688L, 59689L, 59760L, 59761L, 59978L, 59979L, 60650L, 60658L, 60659L, 61278L, 61334L, 61628L, 61629L, 61684L, 61685L, 61753L, 61808L, 61908L, 61964L, 62040L, 62095L, 62162L, 62216L, 62523L, 62524L, 63177L, 63190L, 63195L, 63269L, 63538L, 63594L, 63649L, 63650L, 63875L, 63876L, 63983L, 63984L, 64091L, 64092L, 64200L, 64201L, 64410L, 64465L, 65260L, 65318L, 65379L, 65434L, 65489L, 65544L, 65599L, 65654L, 65709L, 65764L, 65852L, 65982L, 66038L, 66094L, 66150L, 67305L, 67306L, 67414L, 67415L, 67522L, 67523L, 68840L, 69682L, 69737L, 69952L, 69953L, 70010L, 70011L, 70069L, 70070L, 70128L, 70129L, 70681L, 70731L, 70780L, 70805L, 70861L, 70916L, 70971L, 71026L, 71081L, 71136L, 71191L, 71246L, 71721L, 71776L, 71834L, 71892L, 72007L, 72196L, 72943L, 73366L, 73424L, 73482L, 73697L, 73922L, 73948L, 74143L, 74144L, 74154L, 74827L, 75507L, 75687L, 75688L, 75747L, 75748L, 76952L, 76953L, 77012L, 77013L, 77072L, 77073L, 77132L, 77133L, 77192L, 77193L, 77444L, 77445L, 79102L, 79158L, 81759L, 83139L, 83140L, 83190L, 83212L, 83216L, 83234L, 83400L, 83456L, 83512L, 83568L, 83870L, 83871L, 84314L, 84370L, 84881L, 84882L, 84940L, 84941L, 84999L, 85000L, 85179L, 85235L, 86141L, 86142L, 86200L, 86201L, 86259L, 86260L, 86302L, 86328L, 86329L, 86352L, 87205L, 87227L, 87334L, 87335L, 87441L, 87442L, 87550L, 87551L, 88697L, 88755L, 88813L, 89098L, 89546L, 89655L, 89656L, 90115L, 90170L, 90229L, 90231L, 91747L, 91801L, 91855L, 91909L, 91963L, 91964L, 92019L, 92020L, 92075L, 92076L, 92131L, 92132L, 93322L, 93377L, 93432L, 93487L, 94416L, 94472L, 94661L, 94716L, 94834L, 95824L, 95880L, 96510L, 96566L, 96622L, 96678L, 97587L, 97697L, 98518L, 98573L, 98628L, 98683L, 98739L, 98814L, 98822L, 98823L, 98832L, 98833L, 99105L, 99161L, 101081L, 101082L, 101138L, 101139L, 101197L, 101198L, 101254L, 101255L, 101313L, 101369L, 102328L, 102382L, 102437L, 102492L, 102547L, 103503L, 103559L, 104976L, 105011L, 107367L, 107423L, 107479L, 107535L, 107591L, 107646L, 107701L, 107756L, 107894L, 108064L, 108065L, 108121L, 108177L, 108178L, 108235L, 108236L, 108293L, 108294L, 108351L, 108352L, 108409L, 108466L, 108523L, 108524L, 108581L, 108582L, 108639L, 108640L, 108697L, 108698L, 108755L, 108756L, 108878L, 109224L, 109225L, 109283L, 109284L, 109342L, 109343L, 109401L, 109402L, 109460L, 109461L, 109519L, 109520L, 109578L, 109579L, 109637L, 109638L, 109696L, 109697L, 109755L, 109756L, 110154L, 111773L, 111827L, 111881L, 111935L, 113121L, 113176L, 113231L, 113286L, 113305L, 113311L, 113371L, 113377L, 113437L, 113443L, 113503L, 113509L, 113568L, 113574L, 113632L, 113639L, 113697L, 113704L, 113762L, 113769L, 113827L, 113834L, 113892L, 113899L, 113957L, 113964L, 114022L, 114029L, 114087L, 114094L, 114152L, 114159L, 114219L, 114224L, 114349L, 114357L, 114417L, 114425L, 114485L, 114493L, 114553L, 114561L, 114622L, 114630L, 114690L, 114697L, 114757L, 114764L, 114824L, 114831L, 114891L, 114897L, 116021L, 116065L, 116077L, 116122L, 116509L, 116575L, 116637L, 116697L, 116760L, 116822L, 117665L, 117740L, 120115L, 120170L, 120350L, 120351L, 120810L, 120859L, 120922L, 120974L, 121032L, 121088L, 121961L, 121962L, 122073L, 122074L, 122291L, 122292L]]


        # Group PDB ids by species and sequence
        # groups is a list of lists of integer chain identifiers
        self.logger.info("Grouping PDB ids by sequence and species")
        grouper = Grouper(self.config, self.session)
        grouper.use_discrepancy = False
        #grouper.must_enforce_single_species = False
        grouper.must_enforce_single_species = True
        groups = grouper(pdbs)
        if not groups:
            raise core.InvalidState("No groups produced")
        self.logger.info("Grouped PDB ids by sequence and species")


        # only keep chain ids where the chain has at least one nucleotide with a base center and rotation matrix
        # that's a bit silly, we could just wait and look that up later
        has_rotations = ft.partial(self.is_member,
                                   self.known_unit_entries(mod.UnitCenters))
        has_centers = ft.partial(self.is_member,
                                 self.known_unit_entries(mod.UnitRotations))

        # create lists of chain ids from each group, then data method will loop over all pairs
        groups_of_chain_ids = []

        for group in groups:
            chains = it.ifilter(disc.valid_chain, group['members'])
            chains = it.ifilter(has_rotations, chains)
            chains = it.ifilter(has_centers, chains)
            chains = it.imap(op.itemgetter('db_id'), chains)
            # convert chains from iterator to list of chain ids and append to list
            chain_list = sorted(list(chains))
            if len(chain_list) > 1:
                groups_of_chain_ids.append(chain_list)

        # Note how many groups are left now that we filtered as above
        self.logger.info("Found %d groups with at least two chains" % len(groups_of_chain_ids))

        # start with the smallest group
        groups_of_chain_ids.sort(key=len)

        # start with the largest group
        groups_of_chain_ids.sort(key=len,reverse=True)

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


        length1 = len(unit_to_position[chain1])
        length2 = len(unit_to_position[chain2])

        if abs(length1-length2) > 1000 or length1/length2 > 10 or length2/length1 > 10:
            self.logger.warning("Dramatically different number of resolved nucleotides, using discrepancy -1")
            self.logger.info("Chain %s length %d, chain %s length %d" % (chain1,length1,chain2,length2))
        else:

            # Map units in chain2 to their experimental sequence position ids.
            # These are not experimental sequence positions; different sequences
            # have different ids for the same position
            positions2 = []              # list of all positions in chain2
            positions2_to_unit2 = {}     # map those positions back to units
            for unit,position in unit_to_position[chain2].items():
                positions2.append(position)
                positions2_to_unit2[position] = unit
            positions2 = set(positions2) # for faster intersections, I think

            # Loop over units in chain1, map to positions, and intersect with
            # the positions that go with units in chain2
            for unit1,position1 in unit_to_position[chain1].items():
                positions1 = set(position_to_position[position1])
                intersection = positions1 & positions2
                positions2 = positions2 - intersection
                if len(intersection) > 1:
                    self.logger.info("Trouble: Found multiple matches:")
                    for c in intersection:
                        self.logger.info("Matched %s and %s" % (unit1,positions2_to_unit2[c]))
                elif len(intersection) == 1:
                    for c in intersection:
                        unit2 = positions2_to_unit2[c]
                        matching_pairs.append((unit1,unit2))
                else:
                    self.logger.info("No match for %s" % unit1)

            self.logger.info("get_unit_correspondences_intersect: query found %d matching pairs" % len(matching_pairs))

        return matching_pairs


    def get_unit_correspondences(self, corr_id, info1, info2):
        """
        Query the database to find all pairs of unit ids with a given
        correspondence id, matching two given chains and matching PDB ids,
        because that is the fastest query.

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
            exp_map1 = aliased(mod.ExpSeqUnitMapping)           # M1
            exp_map2 = aliased(mod.ExpSeqUnitMapping)           # M2

            mycolumns = [corr_pos.correspondence_id, exp_map1.unit_id.label('unit1'), exp_map2.unit_id.label('unit2')]

            # if mycolumns does not reference the first table to be used in the join, use select_from
            # select_from(corr_pos).\
            # The query without joining units1, units2 takes about 41 seconds on rnatest, for T.th. LSU
            # Joining on units1, units2 to get the pdb id takes about 58 seconds
            # Not filtering for chain takes many, many minutes, so don't do that
            # filtering exp_map1.unit_id.contains(info1['pdb']) takes 46 seconds
            # filtering also on info2['pdb'] takes 27 to 39 seconds, so that seems to be the fastest
            # filtering with like(info1['pdb']+"%") on both 1 and 2 takes 79 seconds

            query = session.query(*mycolumns).\
                join(exp_map1, exp_map1.exp_seq_position_id == corr_pos.exp_seq_position_id_1).\
                join(exp_map2, exp_map2.exp_seq_position_id == corr_pos.exp_seq_position_id_2).\
                filter(corr_pos.correspondence_id == corr_id).\
                filter(exp_map1.chain == info1['chain_name']).\
                filter(exp_map2.chain == info2['chain_name']).\
                filter(exp_map1.unit_id.contains(info1['pdb'])).\
                filter(exp_map2.unit_id.contains(info2['pdb']))

            """
                filter(exp_map1.unit_id.like(info1['pdb']+"%")).\
                filter(exp_map2.unit_id.like(info2['pdb']+"%"))

                limit(10)
                join(units1, units1.unit_id == exp_map1.unit_id).\
                join(units2, units2.unit_id == exp_map2.unit_id).\
                filter(units1.pdb_id == info1['pdb']).\
                filter(units2.pdb_id == info2['pdb'])
                filter(units1.unit_id <> units2.unit_id)
                filter(units1.chain == info1['chain_name']).\
                filter(units2.chain == info2['chain_name'])
            """

            # the query does not seem to actually run until you ask for data from it;
            # setting up is instantaneous, evaluating or counting takes time
            self.logger.info("get_unit_correspondences: query set up for %s" % corr_id)

            matching_pairs = []
            for r in query:
                if r.unit1 and r.unit2 and "|" in r.unit1 and "|" in r.unit2:
                    matching_pairs.append((r.unit1,r.unit2))

            self.logger.info("get_unit_correspondences: query found %d matching pairs" % len(matching_pairs))

            return matching_pairs

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
            # picklefile = 'pickle-FR3D/' + chain_string.replace('|','-') + '_RNA.pickle'
            # chain_chain/comparison.py should read the _NA files instead of the _RNA files, 
            # so they use the glycosidic atom location instead of the base center.  
            # Then hopefully it will still work with RNA and it will also work with DNA.  
            # We will have to check that it looks OK.
            picklefile = 'pickle-FR3D/' + chain_string.replace('|','-') + '_NA.pickle'

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


    def info(self, chain_id):
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
                                  mod.IfeInfo.pdb_id.label('pdb'),
                                  mod.IfeInfo.model,
                                  mod.IfeInfo.ife_id,
                                  ).\
                join(mod.IfeChains,
                     mod.IfeChains.chain_id == mod.ChainInfo.chain_id).\
                join(mod.IfeInfo,
                     mod.IfeInfo.ife_id == mod.IfeChains.ife_id).\
                filter(mod.IfeInfo.new_style == 1).\
                filter(mod.ChainInfo.chain_id == chain_id)

        if not query.count():
            raise core.InvalidState("Could not load chain with id %s" %
                                    chain_id)
        ife = ut.row2dict(query.first())

        with self.session() as session:
            query = session.query(mod.UnitInfo.sym_op,
                                  mod.UnitInfo.alt_id,
                                  ).\
                join(mod.ChainInfo,
                     (mod.ChainInfo.pdb_id == mod.UnitInfo.pdb_id) &
                     (mod.ChainInfo.chain_name == mod.UnitInfo.chain)).\
                join(mod.UnitCenters,
                     mod.UnitCenters.unit_id == mod.UnitInfo.unit_id).\
                join(mod.UnitRotations,
                     mod.UnitRotations.unit_id == mod.UnitInfo.unit_id).\
                filter(mod.ChainInfo.chain_id == chain_id).\
                distinct()

            if not query.count():
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
                filter(info.good_alignment == 1)

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
        seen = set()

        for (unit1,unit2) in unit_pairs:
            if unit1 in seen:
                raise core.InvalidState("gather_matching_centers_rotations: Got duplicate unit1 %s" % unit1)
            seen.add(unit1)

            if unit2 in seen:
                raise core.InvalidState("gather_matching_centers_rotations: Got duplicate unit2 %s" % unit2)
            seen.add(unit2)

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

        if len(c1) < 3:
            self.logger.warning("Too few matched nucleotides to compute discrepancy for %s %s, using magic values instead" %
                                  (info1['name'], info2['name']))
            disc = -1
        else:

            try:
                disc = matrix_discrepancy(c1, r1, c2, r2)

            except Exception as err:
                self.logger.warning("Could not compute discrepancy for %s %s, using magic values instead" %
                                  (info1['name'], info2['name']))
                disc = -1

            if np.isnan(disc):
                self.logger.warning("Could not compute discrepancy for %s %s, using magic values instead" %
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
        self.logger.info("data: Computing discrepancies for a group of %d chains" % L)
        self.logger.info("data: %d discrepancies needed in this group" % (L*(L-1)/2))

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
        self.logger.info("data: Found %d discrepancy values already calculated" % (len(already_computed)/2))

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
        self.logger.info("data: Looking up correspondence ids for %d pairs of chains that need to be computed" % len(check_pairs))

        required_pairs = []
        log_count = 0
        chain_info = {}

        for (chain1_id,chain2_id) in check_pairs:
            corr_id = self.corr_id(chain1_id, chain2_id)
            if corr_id is None:
                if log_count < 20:
                    if not chain1_id in chain_info:
                        chain_info[chain1_id] = self.info(chain1_id)
                    if not chain2_id in chain_info:
                        chain_info[chain2_id] = self.info(chain2_id)
                    info1 = chain_info[chain1_id]
                    info2 = chain_info[chain2_id]

                    self.logger.info("data: No correspondence id between chains %s %s and %s %s" % (info1['ife_id'],chain1_id,info2['ife_id'],chain2_id))
                    log_count = log_count + 1
                else:
                    self.logger.info("data: No correspondence id between chains %s and %s" % (chain1_id,chain2_id))
                # Note: cannot store a discrepancy with a null correspondence id
            else:
                required_pairs.append((corr_id,chain1_id,chain2_id))
        self.logger.info("data: Found %d discrepancy values needing to be calculated" % len(required_pairs))



        # This block is for debugging, to check to see if we get the same discrepancies as before
        Recompute = True   # run the debugging
        Recompute = False

        # check to see that we get the same discrepancy as before for some cases
        # this is for debugging; generally the program will be run with Recompute = False
        if Recompute and len(already_computed_discrepancy) > 0:

            # retrieve chain information for all chains once and store the data
            # 50 seconds for T.th. SSU on rnatest in December 2020

            self.logger.info("data: Retrieving chain information once for each chain")
            chain_info = {}
            for chain_id in chain_ids:
                chain_info[chain_id] = self.info(chain_id)

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

            current = 1
            chain1_seen = set()

            LL = min(len(already_computed_discrepancy),50)

            for (chain1_id,chain2_id,discrepancy) in already_computed_discrepancy[0:LL]:

                # make sure not to repeat any pairings in different order
                if chain1_id < chain2_id:

                    self.logger.info("data: Re-computing discrepancy %d for this group" % (current))
                    current += 1

                    chain1_seen.add(chain1_id)
                    if len(chain1_seen) > 20:
                        chain1_seen = set()
                        allunitdictionary = defaultdict()    # avoid accumulating data forever; reset sometimes

                    info1 = chain_info[chain1_id]
                    info2 = chain_info[chain2_id]

                    # new method
                    self.logger.info("data: Intersect for matching units for chain %s, chain %s" % (info1['ife_id'],info2['ife_id']))
                    new_unit_pairs = self.get_unit_correspondences_intersect(info1,info2,unit_to_position,position_to_position)

                    # filter out units with wrong symmetry or alt id
                    unit_pairs = self.filter_unit_correspondences(new_unit_pairs,info1,info2)

                    # show some matched units to build confidence
                    if len(unit_pairs) > 0:
                        for i in range(0,min(5,len(unit_pairs))):
                            self.logger.info("data: Matched %s and %s" % unit_pairs[i])

                        # load center and rotation data for the current ifes, if not already loaded
                        allunitdictionary = self.load_centers_rotations_pickle(info1,allunitdictionary)
                        allunitdictionary = self.load_centers_rotations_pickle(info2,allunitdictionary)


                    # gather matching centers and rotations for these chains
                    [c1, c2, r1, r2] = self.gather_matching_centers_rotations(unit_pairs,allunitdictionary)
                    self.logger.info("data: Got matching centers and rotations")

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
                self.logger.info("data: Computing discrepancy %d of %d for this group" % (current,len(required_pairs)))

                chain1_seen.add(chain1_id)
                if len(chain1_seen) > 20:
                    chain1_seen = set()
                    allunitdictionary = defaultdict()    # avoid accumulating data forever; reset sometimes

                # Store chain info so that each one is only looked up once.
                # It would seem to be better to just retrieve all of them at once,
                # but that crashed the pipeline once with "QueuePool limit overflow"
                # so it seems better to do them one by one.

                if not chain1_id in chain_info:
                    chain_info[chain1_id] = self.info(chain1_id)
                if not chain2_id in chain_info:
                    chain_info[chain2_id] = self.info(chain2_id)

                info1 = chain_info[chain1_id]
                info2 = chain_info[chain2_id]

                # Future work:
                # Recognize multiple chains for IFEs made of more than one chain; currently only 1st chain is used

                # new method
                self.logger.info("data: Intersect for matching units for chain %s, chain %s" % (info1['ife_id'],info2['ife_id']))
                unit_pairs = self.get_unit_correspondences_intersect(info1,info2,unit_to_position,position_to_position)

                # filter out units with wrong symmetry or alt id
                unit_pairs = self.filter_unit_correspondences(unit_pairs,info1,info2)


                """
                # old method for getting correspondences
                self.logger.info("data: Query for matching units for chain %s, chain %s" % (info1['ife_id'],info2['ife_id']))
                old_unit_pairs = self.get_unit_correspondences(corr_id,info1,info2)

                # filter out units with wrong symmetry or alt id
                old_unit_pairs = self.filter_unit_correspondences(old_unit_pairs,info1,info2)
                self.logger.info("data: Found %d unit id pairs for chain %s chain %s" % (len(unit_pairs),info1['ife_id'],info2['ife_id']))

                self.compare_list_of_pairs(old_unit_pairs,unit_pairs)
                """


                # show some matched units to build confidence
                if len(unit_pairs) > 0:
                    for i in range(0,min(5,len(unit_pairs))):
                        self.logger.info("data: Matched %s and %s" % unit_pairs[i])

                    # load center and rotation data for the current ifes, if not already loaded
                    allunitdictionary = self.load_centers_rotations_pickle(info1,allunitdictionary)
                    allunitdictionary = self.load_centers_rotations_pickle(info2,allunitdictionary)

                # gather matching centers and rotations for these chains
                [c1, c2, r1, r2] = self.gather_matching_centers_rotations(unit_pairs,allunitdictionary)
                self.logger.info("data: Gathered %d matching centers and rotations for %s and %s, %d of %d in this group" % (len(c1),info1['ife_id'],info2['ife_id'],current,len(required_pairs)))

                # compute the discrepancy between these IFEs
                # if wrong numbers of matched nucleotides, discrepancy will be -1
                discrepancies = self.calculate_discrepancy(info1, info2, corr_id, c1, c2, r1, r2)

                self.logger.info("data: Discrepancy to load: %s" % discrepancies[0])
                for d in discrepancies:
                    yield mod.ChainChainSimilarity(**d)
