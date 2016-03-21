"""This is a loader to load the chain to chain discrepancies into the database.
This will look at good correspondences for a given structure and extract all
the aligned chains and then compute the geometric discrepancy between them and
then place them in the database.

This will only compare chains which are an integral part of an IFE.
"""

import itertools as it

import numpy as np
# from sqlalchemy.orm import aliased

from pymotifs import core
from pymotifs import models as mod
import pymotifs.utils as ut
from pymotifs.utils import correspondence as corr
# from pymotifs.utils import temporary_tables as tt

from pymotifs.correspondence.loader import Loader as CorrespondenceLoader
from pymotifs.download import Downloader
from pymotifs.exp_seq.mapping import Loader as ExpSeqUnitMappingLoader
from pymotifs.ife.loader import Loader as IfeLoader

from fr3d.geometry.discrepancy import discrepancy


class MissingAllBaseCenters(Exception):
    """Raised when all bases are missing a base center in some structure.
    """
    pass


class Loader(core.SimpleLoader):
    """A Loader to get all chain to chain similarity data. This will use the
    correspondences between ifes to compute the discrepancy between them.
    """

    allow_no_data = True
    dependencies = set([CorrespondenceLoader, ExpSeqUnitMappingLoader,
                        Downloader, IfeLoader])

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
        return list(pairs)

    def query(self, session, entry, **kwargs):
        """Check if there are any chain_chain_similarity entries for this pdb.
        """

        return session.query(mod.ChainChainSimilarity).\
            filter(mod.ChainChainSimilarity.chain_id_1 == entry[0]).\
            filter(mod.ChainChainSimilarity.chain_id_2 == entry[1]).\
            filter(mod.ChainChainSimilarity.correspondence_id == entry[2])

    def residues(self, cif, model, sym_op, chain, ordering):
        """Get the specified chain and extract only the residues in the
        ordering dictionary and return them in the specified order.

        :cif: The cif object to get the chain from.
        :model: The model number to use.
        :sym_op: The symmetry_operator to use.
        :chain: The name of the chain to access.
        :ordering: The ordering dictonary to use.
        :returns: A sorted list of the residues in the chain.
        """

        residues = cif.residues(model=model, chain=chain, symmetry=sym_op)
        residues = it.ifilter(lambda r: r.unit_id() in ordering, residues)
        residues = sorted(residues, key=lambda r: ordering[r.unit_id()])
        if not residues:
            raise core.InvalidState("Could not find all residues for %s" %
                                    chain)

        return residues

    def discrepancy(self, corr_id, residues1, residues2):
        """Compare the chains. This will filter out all residues in the chain
        that do not have a base center computed.
        """

        valid1 = []
        valid2 = []

        for r1, r2 in zip(residues1, residues2):
            if 'base' not in r1.centers or 'base' not in r2.centers:
                self.logger.debug("Missing base centers for %s, %s", r1, r2)
                continue

            skip = False
            if r1.centers['base'] is None:
                self.logger.warning("Bad center for %s", r1)
                skip = True
            if r2.centers['base'] is None:
                self.logger.warning("Bad center for %s", r2)
                skip = True

            if skip:
                continue

            if len(r1.centers['base']) != len(r2.centers['base']):
                self.logger.warning("Base centers %s, %s differ in size",
                                    r1, r2)
                continue

            valid1.append(r1)
            valid2.append(r2)

        if not valid1:
            self.logger.error("No residues with base centers for first")
            raise MissingAllBaseCenters("Missing first's bases")

        if not valid2:
            self.logger.error("No residues with base centers for second")
            raise MissingAllBaseCenters("Missing second's bases")

        self.logger.debug("Comparing %i pairs of residues", len(residues1))
        disc = discrepancy(valid1, valid2)

        if np.isnan(disc):
            data = (residues1[0].unit_id(), residues2[0].unit_id())
            raise core.InvalidState("NaN for discrepancy using %s, %s" % data)

        return disc

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
                one()
            info = ut.row2dict(result)
            info['model'] = info['model'] or 1

        with self.session() as session:
            info['sym_op'] = session.query(mod.UnitInfo.sym_op).\
                filter_by(pdb_id=info['pdb']).\
                distinct().\
                limit(1).\
                one().\
                sym_op

        return info

    def load(self, ife_id):
        info = self.info(ife_id)
        structure = self.structure(info['pdb'])
        structure.infer_hydrogens()
        return structure, info

    def data(self, entry, **kwargs):
        """Compute all chain to chain similarity data. This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.
        """

        cif1, info1 = self.load(entry[0])
        cif2, info2 = self.load(entry[1])
        corr_id = entry[2]

        util = corr.Helper(self.config, self.session.maker)
        ordering = util.ordering(corr_id, info1['chain_id'], info2['chain_id'])
        if not ordering:
            dat = (corr_id, info1['chain_id'], info2['chain_id'])
            raise core.InvalidState("No ordering for %i, %s, %s", dat)

        residues1 = self.residues(cif1, info1['model'], info1['sym_op'],
                                  info1['chain_name'], ordering)
        residues2 = self.residues(cif2, info2['model'], info2['sym_op'],
                                  info2['chain_name'], ordering)

        if not residues1 or not residues2:
            raise core.Skip("No residues to compare")

        compare = {
            'chain_id_1': info1['chain_id'],
            'chain_id_2': info2['chain_id'],
            'model_1': info1['model'],
            'model_2': info2['model'],
            'correspondence_id': corr_id,
            'discrepancy': self.discrepancy(corr_id, residues1, residues2),
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
