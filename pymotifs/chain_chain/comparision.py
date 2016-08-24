"""This is a loader to load the chain to chain discrepancies into the database.
This will look at good correspondences for a given structure and extract all
the aligned chains and then compute the geometric discrepancy between them and
then place them in the database.

This will only compare chains which are an integral part of an IFE. Each time
it is run it will only add 10 new connections to the database.
"""

import sys
import itertools as it
import functools as ft

import numpy as np
from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod
import pymotifs.utils as ut

from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqUnitMappingLoader
from pymotifs.ife.loader import Loader as IfeLoader
from pymotifs.units.centers import Loader as CenterLoader
from pymotifs.units.rotation import Loader as RotationLoader

from pymotifs.constants import MAX_RESOLUTION_DISCREPANCY
from pymotifs.constants import MIN_NT_DISCREPANCY
from pymotifs.constants import NR_SEQUENCE_ONLY

from fr3d.geometry.discrepancy import matrix_discrepancy


def label_center(table, number):
    return [
        getattr(table, 'unit_id').label('unit%s' % number),
        getattr(table, 'x').label('x%s' % number),
        getattr(table, 'y').label('y%s' % number),
        getattr(table, 'z').label('z%s' % number),
    ]


def label_rotation(table, number):
    return [
        getattr(table, 'cell_0_0').label('cell_00_%s' % number),
        getattr(table, 'cell_0_1').label('cell_01_%s' % number),
        getattr(table, 'cell_0_2').label('cell_02_%s' % number),
        getattr(table, 'cell_1_0').label('cell_10_%s' % number),
        getattr(table, 'cell_1_1').label('cell_11_%s' % number),
        getattr(table, 'cell_1_2').label('cell_12_%s' % number),
        getattr(table, 'cell_2_0').label('cell_20_%s' % number),
        getattr(table, 'cell_2_1').label('cell_21_%s' % number),
        getattr(table, 'cell_2_2').label('cell_22_%s' % number),
    ]


class Loader(core.SimpleLoader):
    """A Loader to get all chain to chain similarity data. This will use the
    correspondences between ifes to compute the discrepancy between them.
    """

    mark = False
    allow_no_data = True
    table = mod.ChainChainSimilarity
    dependencies = set([CorrespondenceLoader, ExpSeqUnitMappingLoader,
                        IfeLoader, CenterLoader, RotationLoader])

    def good_resolution(self, chain):
        return chain['resolution'] is not None and \
            chain['resolution'] <= MAX_RESOLUTION_DISCREPANCY

    def to_process(self, pdbs, **kwargs):
        """This will compute all pairs to compare. This will compute a grouping
        based upon sequence only and then go through

        :params list pdbs: The list of pdbs to use.
        :returns: A list of pairs of chain ids to compare.
        """

        grouper = Grouper(self.config, self.session)
        grouper.use_discrepancy = False
        grouper.must_enforce_single_species = False
        groups = grouper(pdbs)
        if not groups:
            raise core.InvalidState("No groups produced")

        possible = []
        for group in groups:
            chains = it.ifilter(self.good_resolution, group['members'])
            chains = it.imap(op.itemgetter('db_id'), chains)
            possible.extend(it.combinations(chains))
        return possible

    def query(self, session, pair):
        """This will delete *all* entires with the given correspondence id,
        which is likely too many. Better safe than sorry. It's probably ok, I
        think.
        """

        sim = mod.ChainChainSimilarity
        return session.query(sim).\
            filter(((sim.chain_id_1 == pair[0]) & (sim.chain_id_2 == pair[1])) |
                   ((sim.chain_id_1 == pair[1]) & (sim.chain_id_2 == pair[0])))

    def matrices(self, corr_id, info1, info2, name='base'):
        """Load the matrices used to compute discrepancies.

        :param int corr_id: The correspondence id.
        :param dict info1: The result of info for the first chain to compare.
        :param dict info2: The result of info for the second chain to compare.
        :param str name: The type of base center to use.
        :returns: This returns 4 lists, which are in order, centers1, centers2,
        rotation1, rotation2.
        """

        with self.session() as session:
            centers1 = aliased(mod.UnitCenters)
            centers2 = aliased(mod.UnitCenters)
            units1 = aliased(mod.UnitInfo)
            units2 = aliased(mod.UnitInfo)
            rot1 = aliased(mod.UnitRotations)
            rot2 = aliased(mod.UnitRotations)
            corr_units = mod.CorrespondenceUnits

            columns = []
            columns.extend(label_center(centers1, 1))
            columns.extend(label_center(centers2, 2))
            columns.extend(label_rotation(rot1, 1))
            columns.extend(label_rotation(rot2, 2))

            results = session.query(*columns).\
                join(corr_units, corr_units.unit_id_1 == centers1.unit_id).\
                join(centers2, corr_units.unit_id_2 == centers2.unit_id).\
                join(units1, units1.unit_id == corr_units.unit_id_1).\
                join(units2, units2.unit_id == corr_units.unit_id_2).\
                join(rot1, rot1.unit_id == units1.unit_id).\
                join(rot2, rot2.unit_id == units2.unit_id).\
                filter(centers1.name == centers2.name).\
                filter(centers1.name == name).\
                filter(corr_units.correspondence_id == corr_id).\
                filter(units1.pdb_id == info1['pdb']).\
                filter(units1.sym_op == info1['sym_op']).\
                filter(units1.model == info1['model']).\
                filter(units1.chain == info1['chain_name']).\
                filter(units2.pdb_id == info2['pdb']).\
                filter(units2.sym_op == info2['sym_op']).\
                filter(units2.model == info2['model']).\
                filter(units2.chain == info2['chain_name']).\
                filter(rot1.unit_id == centers1.unit_id).\
                filter(rot2.unit_id == centers2.unit_id).\
                filter((units1.alt_id == None) | (units1.alt_id == 'A')).\
                filter((units2.alt_id == None) | (units2.alt_id == 'A')).\
                order_by(corr_units.correspondence_index)

            c1 = []
            c2 = []
            r1 = []
            r2 = []
            seen = set()
            for r in results:
                if r.unit1 in seen:
                    raise core.InvalidState("Got duplicate unit %s" % r.unit1)
                seen.add(r.unit1)

                if r.unit2 in seen:
                    raise core.InvalidState("Got duplicate unit %s" % r.unit2)
                seen.add(r.unit2)

                r1.append(np.array([[r.cell_00_1, r.cell_01_1, r.cell_02_1],
                                    [r.cell_10_1, r.cell_11_1, r.cell_12_1],
                                    [r.cell_20_1, r.cell_21_1, r.cell_22_1]]))
                r2.append(np.array([[r.cell_00_2, r.cell_01_2, r.cell_02_2],
                                    [r.cell_10_2, r.cell_11_2, r.cell_12_2],
                                    [r.cell_20_2, r.cell_21_2, r.cell_22_2]]))

                c1.append(np.array([r.x1, r.y1, r.z1]))
                c2.append(np.array([r.x2, r.y2, r.z2]))

            if not len(c1):
                self.logger.error("No aligned centers for %s", info1['name'])

            if not len(c2):
                self.logger.error("No aligned centers for %s", info2['name'])

            if not len(r1):
                self.logger.error("No aligned rotations for %s", info1['name'])

            if not len(r2):
                self.logger.error("No aligned rotations for %s", info2['name'])

            return np.array(c1), np.array(c2), np.array(r1), np.array(r2)

    def discrepancy(self, corr_id, centers1, centers2, rot1, rot2):
        """Compare the chains. This will filter out all residues in the chain
        that do not have a base center computed.
        """

        self.logger.debug("Comparing %i pairs of residues", len(centers1))
        disc = matrix_discrepancy(centers1, rot1, centers2, rot2)

        if np.isnan(disc):
            raise core.InvalidState("NaN for discrepancy")

        return disc, len(centers1)

    def info(self, chain_id):
        mapping = {}
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
            query = session.query(mod.UnitInfo.sym_op).\
                join(mod.ChainInfo,
                     (mod.ChainInfo.pdb_id == mod.UnitInfo.pdb_id) &
                     (mod.ChainInfo.chain_name == mod.UnitInfo.chain)).\
                filter(mod.ChainInfo.chain_id == chain_id).\
                distinct()

            if not query.count():
                raise core.InvalidState("Could not get sym info for chain %s" %
                                        chain_id)

            possible = set(result.sym_op for result in query)
            for pref in ['1_555', 'P_1']:
                if pref in possible:
                    ife['sym_op'] = pref
                    ife['name'] = ife['ife_id'] + '+' + pref
                    return ife

            sym = sorted(possible)[0]
            ife['sym_op'] = sym
            ife['name'] = ife['ife_id'] + '+' + sym
            return ife

    def __correspondence_query__(self, chain1, chain2):
        with self.session() as session:
            info = mod.CorrespondenceInfo
            mapping1 = aliased(mod.ExpSeqChainMapping)
            mapping2 = aliased(mod.ExpSeqChainMapping)
            query = session.query(info.correspondence_id).\
                join(mapping1, mapping1.exp_seq_id == info.exp_seq_id_1).\
                join(mapping2, mapping2.exp_seq_id == info.exp_seq_id_2).\
                filter(mapping1.chain_id == chain1).\
                filter(mapping2.chain_id == chain2)

            result = query.first()
            if result:
                return result.correspondence_id
            return None

    def corr_id(self, chain_id1, chain_id2):
        """Given two chain ids, load the correspondence id between them. This
        will raise an exception if these chains have not been aligned, or if
        there is more than one alignment between these chains. The chain ids
        should be the ids used in the database.

        :param int chain_id_1: The first chain id.
        :param int chain_id_2: The second chain id.
        :returns: A int for the correspondence id for an alignment between the
        two chains.
        """

        corr_id = self.__correspondence_query__(chain_id1, chain_id2)
        if corr_id is None:
            corr_id = self.__correspondence_query__(chain_id2, chain_id1)

        if corr_id is None:
            raise core.InvalidState("No correspondence for %s %s" %
                                    (chain_id1, chain_id2))
        return corr_id

    def __check_matrices__(self, table, info):
        """Check the that required matrices exist.
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
        """Check if the given info has all the required matrices.
        """

        return self.__check_matrices__(mod.UnitCenters, info) and \
            self.__check_matrices__(mod.UnitRotations, info)

    def entry(self, info1, info2, corr_id):
        """Compute the discrepancy between two given chains.

        :param dict info1: The first chain to use.
        :param dict info2: The second chain to use.
        :param int corr_id: The correspondence id between the chains.
        :returns: A list of dictonaries with the discrepancies betweeen these
        chains.
        """

        if not self.has_matrices(info1):
            self.logger.warning("Missing matrix data for %s", info1['name'])
            return []

        if not self.has_matrices(info2):
            self.logger.warning("Missing matrix data for %s", info2['name'])
            return []

        matrices = self.matrices(corr_id, info1, info2)
        if len(filter(lambda m: len(m), matrices)) != len(matrices):
            self.logger.warning("Did not load all data for %s, %s",
                                info1['name'], info2['name'])
            return []

        try:
            disc, length = self.discrepancy(corr_id, *matrices)
        except Exception as err:
            self.logger.error("Could not compute discrepancy for %s %s" %
                              (info1['name'], info2['name']))
            self.logger.exception(err)
            return []

        compare = {
            'chain_id_1': info1['chain_id'],
            'chain_id_2': info2['chain_id'],
            'model_1': info1['model'],
            'model_2': info2['model'],
            'correspondence_id': corr_id,
            'discrepancy': disc,
            'num_nucleotides': length,
        }

        reversed = dict(compare)
        reversed['chain_id_1'] = compare['chain_id_2']
        reversed['chain_id_2'] = compare['chain_id_1']
        reversed['model_1'] = compare['model_2']
        reversed['model_2'] = compare['model_1']

        return [
            compare,
            reversed
        ]

    def data(self, pair, **kwargs):
        """Compute all chain to chain similarity data. This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.
        """

        chain1, chain2 = pair
        info1 = self.info(chain1)
        info2 = self.info(chain2)
        corr_id = self.corr_id(chain1, chain2)
        return self.entry(info1, info2, corr_id)
