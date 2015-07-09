"""This is a loader to load the chain to chain comparisons into the database.
This will look at good correspondences for a given structure and extract all
the aligned chains and then compute the geometric discrepancy between them and
then place them in the database.
"""

import itertools as it

import numpy as np

from pymotifs import core
from pymotifs.utils import correspondence as corr

from pymotifs.models import ChainInfo
from pymotifs.models import ChainChainSimilarity as Similarity
from pymotifs.correspondence import Loader as CorrespondenceLoader
from pymotifs.download import Downloader
from pymotifs.exp_seq.mapping import Loader as ExpSeqUnitMappingLoader
from pymotifs.autonomous import Loader as AutonomousLoader

from fr3d.geometry.discrepancy import discrepancy


class MissingAllBaseCenters(Exception):
    pass


class MissingReference(Exception):
    pass


class Loader(core.Loader):
    """A Loader to get all chain to chain similarity data.
    """

    allow_no_data = True
    dependencies = set([CorrespondenceLoader, ExpSeqUnitMappingLoader,
                        Downloader, AutonomousLoader])

    def known(self, pdb):
        with self.session() as session:
            query = session.query(Similarity.correspondence_id,
                                  Similarity.chain_id1,
                                  Similarity.chain_id2,
                                  ChainInfo.pdb_id).\
                join(ChainInfo, ChainInfo.id == Similarity.chain_id1).\
                filter(ChainInfo.pdb_id == pdb)

            known = []
            for result in query:
                known.append((result.correspondence_id, result.chain_id1,
                              result.pdb_id, result.chain_id2))
        return known

    def possible(self, pdb):
        data = []
        util = corr.Helper(self.config, self.session.maker)
        others = util.pdbs(pdb)
        if not others:
            raise core.Skip("No pdbs to compare %s to" % pdb)

        for pdb2 in others:

            corresponding = util.chains(pdb, pdb2)
            if not corresponding:
                self.logger.error("Could not find chains to compare in %s %s",
                                  pdb, pdb2)

            for corr_id, chain1, chain2 in corresponding:
                data.append((corr_id, chain1['id'], chain2['pdb'],
                             chain2['id']))
        return data

    def missing(self, pdb):
        return set(self.possible(pdb)) - set(self.known(pdb))

    def has_data(self, pdb, **kwargs):
        """Check if there are any chain_chain_similarity entries for this pdb.
        """
        return not bool(self.missing(pdb))

    def remove(self, pdb, **kwargs):
        """Remove all chain_chain_similarity entries for this pdb.
        """

        with self.session() as session:
            query = session.query(Similarity).\
                join(ChainInfo, ChainInfo.id == Similarity.chain_id1).\
                filter(ChainInfo.pdb_id == pdb)

            ids = [result.id for result in query]

        if not ids:
            self.logger.info("Nothing to remove for %s", pdb)
            return None

        with self.session() as session:
            session.query(Similarity).\
                filter(Similarity.id.in_(ids)).\
                delete(synchronize_session=False)

    def residues(self, cif, name, ordering):
        """Get the specified chain and extract only the residues in the
        ordering dictionary and return them in the specified order.

        :cif: The cif object to get the chain from.
        :chain: The name of the chain to access.
        :returns: A sorted list of the residues in the chain.
        """

        residues = cif.residues(model=1, chain=name)
        residues = it.ifilter(lambda r: r.unit_id() in ordering, residues)
        residues = sorted(residues, key=lambda r: ordering[r.unit_id()])
        if not residues:
            self.logger.error("Could not find any residues for %s", name)
            return None

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
            raise MissingAllBaseCenters()

        if not valid2:
            self.logger.error("No residues with base centers for second")
            disc = None

        self.logger.debug("Comparing %i pairs of residues", len(residues1))

        try:
            disc = discrepancy(valid1, valid2)
        except Exception as err:
            self.logger.error("Error computing discrepancy")
            self.logger.exception(err)
            disc = None

        if np.isnan(disc):
            self.logger.error("Got NaN for discrepancy between")
            return None

        return disc

    def chain_name(self, chain_id):
        with self.session() as session:
            return session.query(ChainInfo).\
                filter_by(id=chain_id).\
                one().chain_name

    def data(self, pdb, **kwargs):
        """Compute all chain to chain similarity data. This will get all
        corresponding chains to chain alignment for all chains in this pdb and
        determine the geometric similarity for the chains, if the alignment is
        a good one.
        """

        data = []
        cif1 = self.structure(pdb)
        cif1.infer_hydrogens()
        missing = self.missing(cif1.pdb)
        grouped = it.groupby(missing, lambda m: m[2])
        util = corr.Helper(self.config, self.session.maker)

        for index, (pdb2, corresponding) in enumerate(grouped):
            self.logger.debug("Comparing to %s %i", pdb2, index)

            cif2 = self.structure(pdb2)
            cif2.infer_hydrogens()

            corresponding = list(corresponding)
            for corr_id, chain_id1, pdb2, chain_id2 in corresponding:
                self.logger.debug("Comparing %s to %s", chain_id1, chain_id2)

                ordering = util.ordering(corr_id, chain_id1, chain_id2)
                if not ordering:
                    self.logger.error("No ordering for %i, %s, %s",
                                      corr_id, chain_id1, chain_id2)
                    continue

                chain_name1 = self.chain_name(chain_id1)
                chain_name2 = self.chain_name(chain_id2)
                residues1 = self.residues(cif1, chain_name1, ordering)
                residues2 = self.residues(cif2, chain_name2, ordering)

                if not residues1 or not residues2:
                    self.logger.error("No residues to compare")
                    continue

                try:
                    compare = {
                        'chain_id1': chain_id1,
                        'chain_id2': chain_id2,
                        'correspondence_id': corr_id,
                        'discrepancy': self.discrepancy(corr_id, residues1,
                                                        residues2)
                    }
                except:
                    continue

                reversed = dict(compare)
                reversed['chain_id1'] = compare['chain_id2']
                reversed['chain_id2'] = compare['chain_id1']
                data.append(Similarity(**compare))
                data.append(Similarity(**reversed))

        return data
