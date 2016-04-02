"""This is a loader to load the chain to chain discrepancies into the database.
This will look at good correspondences for a given structure and extract all
the aligned chains and then compute the geometric discrepancy between them and
then place them in the database.

This will only compare chains which are an integral part of an IFE.
"""

import itertools as it
import collections as coll

import numpy as np
from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod
import pymotifs.utils as ut
from pymotifs.constants import NR_DISCREPANCY_CUTOFF
from pymotifs.utils import correspondence as corr

from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.exp_seq.mapping import Loader as ExpSeqUnitMappingLoader
from pymotifs.ife.loader import Loader as IfeLoader
from pymotifs.units.centers import Loader as CenterLoader
from pymotifs.units.rotation import Loader as RotationLoader

from fr3d.geometry.discrepancy import matrix_discrepancy


class Loader(core.SimpleLoader):
    """A Loader to get all chain to chain similarity data. This will use the
    correspondences between ifes to compute the discrepancy between them.
    """

    mark = False
    dependencies = set([CorrespondenceLoader, ExpSeqUnitMappingLoader,
                        IfeLoader, CenterLoader, RotationLoader])
    max_new_connections = 10

    def __init__(self, *args, **kwargs):
        super(Loader, self).__init__(*args, **kwargs)
        self.new_updates = coll.defaultdict(int)

    def ifes_info(self, pdb):
        with self.session() as session:
            query = session.query(mod.IfeInfo.ife_id,
                                  mod.IfeChains.chain_id).\
                join(mod.IfeChains,
                     mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
                filter(mod.IfeInfo.pdb_id == pdb).\
                filter(mod.IfeChains.index == 0)
            return [ut.row2dict(result) for result in query]

    def find_correspondence(self, pair):
        with self.session() as session:
            pdbs = mod.CorrespondencePdbs
            info = mod.CorrespondenceInfo
            result = session.query(pdbs.correspondence_id,
                                   info.good_alignment).\
                join(info,
                     info.correspondence_id == pdbs.correspondence_id).\
                filter(pdbs.chain_id_1 == pair[0]['chain_id']).\
                filter(pdbs.chain_id_2 == pair[1]['chain_id']).\
                first()

            return (pair[0], pair[1], ut.row2dict(result))

    def to_process(self, pdbs, **kwargs):
        """This will figure out all chains that need to be processed to
        determine the discrepancies.
        """

        def as_tuple(pair):
            return (pair[0]['ife_id'], pair[1]['ife_id'],
                    pair[2]['correspondence_id'])

        ifes = it.imap(self.ifes_info, pdbs)
        ifes = it.chain.from_iterable(ifes)
        pairs = it.combinations(ifes, 2)
        pairs = it.imap(self.find_correspondence, pairs)
        pairs = it.ifilter(lambda p: p[2], pairs)
        pairs = it.ifilter(lambda p: p[2]['good_alignment'], pairs)
        pairs = it.imap(as_tuple, pairs)
        return sorted(pairs, key=lambda p: (p[0], p[1]))

    def has_data(self, entry, **kwargs):
        if super(Loader, self).has_data(entry, **kwargs):
            return True
        return self.new_updates[entry[0]] > self.max_new_connections

    def query(self, session, entry, **kwargs):
        """Check if there are any chain_chain_similarity entries for this pdb.
        """

        ife1 = aliased(mod.IfeChains)
        ife2 = aliased(mod.IfeChains)
        cc = mod.ChainChainSimilarity
        return session.query(cc).\
            join(ife1, ife1.chain_id == cc.chain_id_1).\
            join(ife2, ife2.chain_id == cc.chain_id_2).\
            filter(ife1.ife_id == entry[0]).\
            filter(ife2.ife_id == entry[1]).\
            filter(ife1.index == 0).\
            filter(ife2.index == 0).\
            filter(cc.correspondence_id == entry[2])

    def centers(self, corr_id, info1, info2, ordering, name='base'):
        with self.session() as session:
            centers1 = aliased(mod.UnitCenters)
            centers2 = aliased(mod.UnitCenters)
            units1 = aliased(mod.UnitInfo)
            units2 = aliased(mod.UnitInfo)
            corr = mod.CorrespondenceUnits
            results = session.query(centers1.unit_id.label('unit1'),
                                    centers1.x.label('x1'),
                                    centers1.y.label('y1'),
                                    centers1.z.label('z1'),
                                    centers2.unit_id.label('unit2'),
                                    centers2.x.label('x2'),
                                    centers2.y.label('y2'),
                                    centers2.z.label('z2')).\
                join(corr, corr.unit_id_1 == centers1.unit_id).\
                join(centers2, corr.unit_id_2 == centers2.unit_id).\
                join(units1, units1.unit_id == corr.unit_id_1).\
                join(units2, units2.unit_id == corr.unit_id_2).\
                filter(centers1.name == centers2.name).\
                filter(centers1.name == name).\
                filter(corr.correspondence_id == corr_id).\
                filter(units1.pdb_id == info1['pdb']).\
                filter(units1.sym_op == info1['sym_op']).\
                filter(units1.model == info1['model']).\
                filter(units1.chain == info1['chain_name']).\
                filter(units2.pdb_id == info2['pdb']).\
                filter(units2.sym_op == info2['sym_op']).\
                filter(units2.model == info2['model']).\
                filter(units2.chain == info2['chain_name']).\
                order_by(corr.correspondence_index).\
                all()

            results.sort(key=lambda r: ordering[r.unit1])

            first = []
            second = []
            seen = set()
            for result in results:
                if result.unit1 in seen:
                    raise core.InvalidState("Got duplicate unit %s" % unit1)
                if result.unit2 in seen:
                    raise core.InvalidState("Got duplicate unit %s" % unit2)

                seen.add(result.unit1)
                seen.add(result.unit2)

                first.append(np.array([result.x1, result.y1, result.z1]))
                second.append(np.array([result.x2, result.y2, result.z2]))
            return first, second

    def rotations(self, corr_id, info1, info2, ordering):
        with self.session() as session:
            rot1 = aliased(mod.UnitRotations)
            rot2 = aliased(mod.UnitRotations)
            units1 = aliased(mod.UnitInfo)
            units2 = aliased(mod.UnitInfo)
            corr = mod.CorrespondenceUnits
            results = session.query(rot1.unit_id.label('unit1'),
                                    rot1.cell_0_0.label('f_cell_0_0'),
                                    rot1.cell_0_1.label('f_cell_0_1'),
                                    rot1.cell_0_2.label('f_cell_0_2'),
                                    rot1.cell_1_0.label('f_cell_1_0'),
                                    rot1.cell_1_1.label('f_cell_1_1'),
                                    rot1.cell_1_2.label('f_cell_1_2'),
                                    rot1.cell_2_0.label('f_cell_2_0'),
                                    rot1.cell_2_1.label('f_cell_2_1'),
                                    rot1.cell_2_2.label('f_cell_2_2'),
                                    rot2.unit_id.label('unit2'),
                                    rot2.cell_0_0.label('s_cell_0_0'),
                                    rot2.cell_0_1.label('s_cell_0_1'),
                                    rot2.cell_0_2.label('s_cell_0_2'),
                                    rot2.cell_1_0.label('s_cell_1_0'),
                                    rot2.cell_1_1.label('s_cell_1_1'),
                                    rot2.cell_1_2.label('s_cell_1_2'),
                                    rot2.cell_2_0.label('s_cell_2_0'),
                                    rot2.cell_2_1.label('s_cell_2_1'),
                                    rot2.cell_2_2.label('s_cell_2_2')).\
                join(corr, corr.unit_id_1 == rot1.unit_id).\
                join(rot2, corr.unit_id_2 == rot2.unit_id).\
                join(units1, rot1.unit_id == units1.unit_id).\
                join(units2, rot2.unit_id == units2.unit_id).\
                filter(units1.pdb_id == info1['pdb']).\
                filter(units1.sym_op == info1['sym_op']).\
                filter(units1.model == info1['model']).\
                filter(units1.chain == info1['chain_name']).\
                filter(units2.pdb_id == info2['pdb']).\
                filter(units2.sym_op == info2['sym_op']).\
                filter(units2.model == info2['model']).\
                filter(units2.chain == info2['chain_name']).\
                filter(corr.correspondence_id == corr_id).\
                order_by(corr.correspondence_index).\
                all()

            # results.sort(key=lambda r: ordering[r.unit1])

            r1 = []
            r2 = []
            seen = set()
            for r in results:
                if r.unit1 in seen:
                    raise core.InvalidState("Got duplicate unit %s" % r.unit1)
                if r.unit2 in seen:
                    raise core.InvalidState("Got duplicate unit %s" % r.unit2)

                seen.add(r.unit1)
                seen.add(r.unit2)

                r1.append(np.array([[r.f_cell_0_0, r.f_cell_0_1, r.f_cell_0_2],
                                    [r.f_cell_1_0, r.f_cell_1_1, r.f_cell_1_2],
                                    [r.f_cell_2_0, r.f_cell_2_1, r.f_cell_2_2]]))
                r2.append(np.array([[r.s_cell_0_0, r.s_cell_0_1, r.s_cell_0_2],
                                    [r.s_cell_1_0, r.s_cell_1_1, r.s_cell_1_2],
                                    [r.s_cell_2_0, r.s_cell_2_1, r.s_cell_2_2]]))
            return r1, r2


    def discrepancy(self, corr_id, centers1, centers2, rot1, rot2):
        """Compare the chains. This will filter out all residues in the chain
        that do not have a base center computed.
        """

        if not centers1 or not rot1:
            self.logger.error("No residues with base centers/rot for first")
            return -1, 0

        if not centers2 or not rot2:
            self.logger.error("No residues with base centers/rot for second")
            return -1, 0

        self.logger.debug("Comparing %i pairs of residues", len(centers1))
        disc = matrix_discrepancy(centers1, centers2, rot1, rot2)

        if np.isnan(disc):
            raise core.InvalidState("NaN for discrepancy")

        return disc, len(centers1)

    def info(self, ife_id):
        with self.session() as session:
            result = session.query(mod.IfeInfo.pdb_id.label('pdb'),
                                   mod.ChainInfo.chain_name,
                                   mod.ChainInfo.chain_id,
                                   mod.IfeInfo.model).\
                join(mod.IfeChains,
                     mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
                join(mod.ChainInfo,
                     mod.ChainInfo.chain_id == mod.IfeChains.chain_id).\
                filter(mod.IfeInfo.ife_id == ife_id).\
                filter(mod.IfeInfo.new_style == 1).\
                filter(mod.IfeChains.index == 0).\
                one()
            info = ut.row2dict(result)

        with self.session() as session:
            info['sym_op'] = session.query(mod.UnitInfo.sym_op).\
                filter_by(pdb_id=info['pdb']).\
                distinct().\
                limit(1).\
                one().\
                sym_op

        return info

    def ordering(self, corr_id, info1, info2):
        util = corr.Helper(self.config, self.session.maker)
        ordering = util.ordering(corr_id, info1['chain_id'], info2['chain_id'])
        if not ordering:
            dat = (corr_id, info1['chain_id'], info2['chain_id'])
            raise core.InvalidState("No ordering for %i, %s, %s", dat)
        return ordering

    def data(self, entry, **kwargs):
        """Compute all chain to chain similarity data. This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.
        """

        info1 = self.info(entry[0])
        info2 = self.info(entry[1])
        corr_id = entry[2]

        ordering = self.ordering(corr_id, info1, info2)
        centers1, centers2 = self.centers(corr_id, info1, info2, ordering)
        rot1, rot2 = self.rotations(corr_id, info1, info2, ordering)
        disc, length = self.discrepancy(corr_id, centers1, centers2, rot1, rot2)

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

        if discrepancy <= NR_DISCREPANCY_CUTOFF:
            self.new_updates[entry[0]] += 1

        return [
            mod.ChainChainSimilarity(**compare),
            mod.ChainChainSimilarity(**reversed)
        ]
