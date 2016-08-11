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
from pymotifs.nr.builder import Loader as ReleaseLoader

from pymotifs.constants import MAX_RESOLUTION_DISCREPANCY
from pymotifs.constants import MIN_NT_DISCREPANCY

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


class Loader(core.Loader):
    """A Loader to get all chain to chain similarity data. This will use the
    correspondences between ifes to compute the discrepancy between them.
    """

    mark = False
    allow_no_data = True
    dependencies = set([CorrespondenceLoader, ExpSeqUnitMappingLoader,
                        ReleaseLoader, IfeLoader, CenterLoader,
                        RotationLoader])

    def to_process(self, pdbs, **kwargs):
        with self.session() as session:
            info = mod.CorrespondenceInfo
            query = session.query(info).\
                filter(info.good_alignment == True)
            return [result.correspondence_id for result in query]

    def remove(self, corr_id, **kwargs):
        """This will delete *all* entires with the given correspondence id,
        which is likely too many. Better safe than sorry. It's probably ok, I
        think.
        """

        with self.session() as session:
            session.query(mod.ChainChainSimilarity).\
                filter_by(correspondence_id=corr_id).\
                delete(synchronize_session=False)

    def has_data(self, *args, **kwargs):
        """We always pretend as if we have no data for this correspondence id
        as it is can take a long time to compute if we have any chains that
        need to be aligned. Thus we pretend we always need to compute stuff and
        then allow this stage to produce no data.
        """
        return False

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

    def ifes(self, exp_seq_id):
        mapping = {}
        with self.session() as session:
            exp_mapping = mod.ExpSeqChainMapping
            query = session.query(mod.IfeInfo.pdb_id.label('pdb'),
                                  mod.ChainInfo.chain_name,
                                  mod.ChainInfo.chain_id,
                                  mod.IfeInfo.model,
                                  mod.IfeInfo.ife_id,
                                  ).\
                join(mod.IfeChains,
                     mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
                join(mod.ChainInfo,
                     mod.ChainInfo.chain_id == mod.IfeChains.chain_id).\
                join(exp_mapping,
                     exp_mapping.chain_id == mod.IfeChains.chain_id).\
                join(mod.PdbInfo,
                     mod.PdbInfo.pdb_id == mod.IfeInfo.pdb_id).\
                filter(mod.IfeInfo.new_style == 1).\
                filter(mod.PdbInfo.resolution != None).\
                filter(mod.PdbInfo.resolution <= MAX_RESOLUTION_DISCREPANCY).\
                filter(mod.ExpSeqChainMapping.exp_seq_id == exp_seq_id)

            pdbs = set()
            for result in query:
                pdbs.add(result.pdb)
                mapping[result.chain_id] = ut.row2dict(result)

        if not mapping:
            self.logger.warning("Loaded no possible ifes in: %s", exp_seq_id_1)
            return []

        with self.session() as session:
            query = session.query(mod.UnitInfo.sym_op,
                                  mod.ChainInfo.chain_id).\
                join(mod.ChainInfo,
                     (mod.ChainInfo.pdb_id == mod.UnitInfo.pdb_id) &
                     (mod.ChainInfo.chain_name == mod.UnitInfo.chain)).\
                filter(mod.UnitInfo.pdb_id.in_(pdbs)).\
                distinct()

            for result in query:
                if result.chain_id in mapping:
                    current = mapping[result.chain_id]
                    current['sym_op'] = result.sym_op
                    current['name'] = current['ife_id'] + '+' + result.sym_op

        return mapping.values()

    def known(self, corr_id):
        with self.session() as session:
            query = session.query(mod.ChainChainSimilarity).\
                filter_by(correspondence_id=corr_id)
            return set([(r.chain_id_1, r.chain_id_2) for r in query])

    def is_known(self, known, pair):
        key = (pair[0]['chain_id'], pair[1]['chain_id'])
        return key in known

    def aligned_ifes(self, corr_id):
        with self.session() as session:
            corr = session.query(mod.CorrespondenceInfo).\
                filter_by(correspondence_id=corr_id).\
                one()
            exp1, exp2 = corr.exp_seq_id_1, corr.exp_seq_id_2
        ife1 = self.ifes(exp1)
        ife2 = self.ifes(exp2)

        possible = it.product(ife1, ife2)
        known = self.known(corr_id)
        is_known = ft.partial(self.is_known, known)
        pairs = it.ifilterfalse(is_known, possible)
        pairs = it.ifilter(lambda p: p[0] != p[1], pairs)
        pairs = list(pairs)

        if not len(pairs):
            raise core.Skip("Done all chains")

        self.logger.info("Found %i possible ife pairs to compare", len(pairs))

        return pairs

    def __check_matrices__(self, table, info):
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
        return self.__check_matrices__(mod.UnitCenters, info) and \
            self.__check_matrices__(mod.UnitRotations, info)

    def entry(self, info1, info2, corr_id):
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
            mod.ChainChainSimilarity(**compare),
            mod.ChainChainSimilarity(**reversed)
        ]

    def data(self, corr_id, **kwargs):
        """Compute all chain to chain similarity data. This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.
        """

        data = []
        ifes = self.aligned_ifes(corr_id)
        for index, (ife1, ife2) in enumerate(ifes):
            self.logger.info("Comparing %i: %s, %s",
                             index, ife1['name'], ife2['name'])
            data.extend(self.entry(ife1, ife2, corr_id))
        return data
