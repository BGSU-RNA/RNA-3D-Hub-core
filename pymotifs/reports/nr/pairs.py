"""This module contains the logic to create the NR reports.
"""

import operator as op

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict
from pymotifs.utils.discrepancy import should_compare_chain_discrepancy
from pymotifs.constants import MIN_NT_DISCREPANCY
from pymotifs.constants import MAX_RESOLUTION_DISCREPANCY

from sqlalchemy.orm import aliased


class Pairs(core.Reporter):
    headers = [
        'Group',
        'Release',
        'IFE1',
        'IFE2',
        'Discrepancy',
        'Alignment',
        'Good Alignment',
    ]

    def data_query(self, release, resolution):
        classes = mod.NrClasses
        chains1 = aliased(mod.NrChains)
        chains2 = aliased(mod.NrChains)
        ife1 = aliased(mod.IfeChains)
        ife2 = aliased(mod.IfeChains)
        mapping1 = aliased(mod.ExpSeqChainMapping)
        mapping2 = aliased(mod.ExpSeqChainMapping)
        ccs = mod.ChainChainSimilarity
        corr = aliased(mod.CorrespondenceInfo)
        rev_corr = aliased(mod.CorrespondenceInfo)
        exp1 = aliased(mod.ExpSeqInfo)
        exp2 = aliased(mod.ExpSeqInfo)
        pdb1 = aliased(mod.PdbInfo)
        pdb2 = aliased(mod.PdbInfo)
        chain_info1 = aliased(mod.ChainInfo)
        chain_info2 = aliased(mod.ChainInfo)
        with self.session() as session:
            query = session.query(classes.name.label('Group'),
                                  classes.nr_release_id.label('Release'),
                                  chains1.ife_id.label('IFE1'),
                                  chains2.ife_id.label('IFE2'),
                                  ccs.discrepancy.label('Discrepancy'),
                                  corr.match_count.label('forward_match'),
                                  rev_corr.match_count.label('rev_match'),
                                  corr.good_alignment.label('forward_good'),
                                  rev_corr.good_alignment.label('rev_good'),
                                  exp1.length.label('len1'),
                                  exp2.length.label('len2'),
                                  pdb1.resolution.label('res1'),
                                  pdb2.resolution.label('res2'),
                                  pdb1.experimental_technique.label('meth1'),
                                  pdb2.experimental_technique.label('meth2'),
                                  ).\
                join(chains1, chains1.nr_class_id == classes.nr_class_id).\
                join(chains2, chains2.nr_class_id == classes.nr_class_id).\
                join(ife1,
                     (ife1.ife_id == chains1.ife_id) & (ife1.index == 0)
                     ).\
                join(ife2,
                     (ife2.ife_id == chains2.ife_id) & (ife2.index == 0)
                     ).\
                join(chain_info1, chain_info1.chain_id == ife1.chain_id).\
                join(chain_info2, chain_info2.chain_id == ife2.chain_id).\
                join(pdb1, pdb1.pdb_id == chain_info1.pdb_id).\
                join(pdb2, pdb2.pdb_id == chain_info2.pdb_id).\
                join(mapping1, mapping1.chain_id == ife1.chain_id).\
                join(mapping2, mapping2.chain_id == ife2.chain_id).\
                join(exp1, exp1.exp_seq_id == mapping1.exp_seq_id).\
                join(exp2, exp2.exp_seq_id == mapping2.exp_seq_id).\
                outerjoin(ccs,
                          (ccs.chain_id_1 == ife1.chain_id) &
                          (ccs.chain_id_2 == ife2.chain_id)
                          ).\
                outerjoin(corr,
                          (corr.exp_seq_id_1 == mapping1.exp_seq_id) &
                          (corr.exp_seq_id_2 == mapping2.exp_seq_id)
                          ).\
                outerjoin(rev_corr,
                          (rev_corr.exp_seq_id_1 == mapping2.exp_seq_id) &
                          (rev_corr.exp_seq_id_2 == mapping1.exp_seq_id)
                          ).\
                filter(chains1.nr_chain_id != chains2.nr_chain_id).\
                filter(ife1.chain_id != ife2.chain_id).\
                filter(classes.resolution == resolution).\
                filter(classes.nr_release_id == release).\
                distinct().\
                order_by(classes.name, ife1.ife_id, ife2.ife_id)
            return [row2dict(e) for e in query]

    def as_chain(self, n, entry):
        return {
            'length': entry['len%s' % n],
            'method': entry['meth%s' % n],
            'resolution': entry['res%s' % n]
        }

    def update_discrepancy(self, entry):
        chain1 = self.as_chain(1, entry)
        chain2 = self.as_chain(2, entry)
        if not should_compare_chain_discrepancy(chain1) or \
                not should_compare_chain_discrepancy(chain2):
            entry['Discrepancy'] = None

    def update_alignment(self, entry):
        for_match = entry['forward_match']
        rev_match = entry['rev_match']
        len1 = entry['len1']
        len2 = entry['len2']
        for_good = entry['forward_good']
        rev_good = entry['rev_good']
        entry['Alignment'] = None
        entry['Good Alignment'] = int(for_good or rev_good)

        if for_match is not None or rev_match is not None:
            match = float(for_match or rev_match)
            entry['Alignment'] = match / float(min(len1, len2))

    def data(self, entries, **kwargs):
        """Create a report about pairs in the NR set.

        Parameters
        ----------
        release : str
            The NR release id to use
        resolution : str
            Resolution cutoff for groups to use.

        Yields
        --------
        entry : dict
            A dict that can written for the report.
        """

        for entry in self.data_query(*entries):
            self.update_discrepancy(entry)
            self.update_alignment(entry)
            del entry['res1']
            del entry['res2']
            del entry['forward_match']
            del entry['rev_match']
            del entry['len1']
            del entry['len2']
            del entry['meth1']
            del entry['meth2']
            del entry['forward_good']
            del entry['rev_good']
            yield entry
