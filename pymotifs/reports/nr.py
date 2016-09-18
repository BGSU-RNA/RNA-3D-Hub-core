"""This module contains the logic to create the NR reports.
"""

import operator as op
import itertools as it
import collections as coll

from pymotifs import models as mod
from pymotifs.utils import row2dict

from sqlalchemy.orm import aliased
from sqlalchemy.sql.expression import func

Report = coll.namedtuple('Report', ['headers', 'rows'])


"""The list of headers to write in the report."""
GROUP_HEADERS = [
    'Group',
    'Rank',
    'Type',
    'IFE id',
    'PDB',
    'Chains',
    'Protein Species',
    'Protein Compound',
    'RNA Species',
    'Nucleic Acid Compound',
    'Name',
    'Observed Length',
    'Experimental Length',
    'Experimental Sequence',
]

PAIRS_HEADERS = [
    'Group',
    'Release',
    'IFE1',
    'IFE2',
    'Discrepancy',
    'Alignment',
]

compound = op.itemgetter('compound')
species = op.itemgetter('source')
macro_type = op.itemgetter('entity_macromolecule_type')
is_type = lambda t: lambda c: macro_type(c) == t
is_protein = is_type('Polypeptide(L)')
is_rna = is_type('Polyribonucleotide (RNA)')

PdbInfo = coll.namedtuple('PdbInfo', ['protein_chains', 'rna_chains'])
ChainInfo = coll.namedtuple('ChainInfo', ['name', 'species', 'compound',
                                          'observed', 'experimental',
                                          'sequence'])


class Entry(object):
    """A class to represent the

    Attributes
    ----------
    group : str
        The name of the equivlance class this
    rank : int
        The rank of the IFE in it's class
    ife_id : str
        The IFE id.
    pdb_id : str
        The PDB id of the ife.
    chains : list
        A list of chain ids, that are a part of this IFE id.
    names : list
        List of chain names in this IFE.
    proteins : list
        List of protein names in this structure.
    protein_species : list
        List of species assignments for proteins in this structure.
    species : list
        List of species assignments for the RNA chains in this IFE.
    compound : str
        Name of the largest chain in the IFE.
    observed : int
        Number of observed residues in the ife
    experimental : int
        Number of residues in the experimental sequence
    sequence : str
        The experimental sequence
    """

    def __init__(self, group, rank, ife_id, pdb_id, chain_ids):
        """Create a new `Entity`.

        Parameters
        ----------
        group : str
            The name of the equivlance class this
        rank : int
            The rank of the IFE in it's class
        ife_id : str
            The IFE id.
        pdb_id : str
            The PDB id of the ife.
        chain_ids : list
            A list of chain ids, that are a part of this IFE id.
        """

        self.group = group
        self.rank = rank
        self.ife_id = ife_id
        self.pdb_id = pdb_id
        self.chains = chain_ids
        self.names = []
        self.proteins = []
        self.protein_species = []
        self.species = None
        self.compound = ''
        self.observed = 0
        self.experimental = 0
        self.sequence = ''

    def add_info(self, info):
        """Use the given `Info` object to update the attributes of this object.
        This will update `proteins`, `protein_species`, `species`, `compound`,
        `observed`, `experimental` and `sequence` attributes.

        Parameters
        ----------
        info : Info
            The info object to use.
        """

        pdb_info = info.pdbs[self.pdb_id]
        prot_chains = pdb_info.protein_chains
        self.proteins = [info.chains[c].compound for c in prot_chains]
        self.protein_species = [info.chains[c].species for c in prot_chains]

        for index, chain_id in enumerate(self.chains):
            chain_info = info.chains[chain_id]
            self.names.append(chain_info.name)
            if index == 0:
                self.species = chain_info.species
                self.compound = chain_info.compound
                self.observed = chain_info.observed
                self.experimental = chain_info.experimental
                self.sequence = chain_info.sequence

    @property
    def type(self):
        if self.rank == 0:
            return 'rep'
        return 'member'

    def named(self, max_proteins=5):
        """Turn this `Entity` into a dictonary with named columns for the
        report.

        Parameters
        ----------
        max_proteins : 5, optional
            The maximal number of proteins in the structure to allow before not
            writing protein level information. If set to None, this will always
            write protein information.

        Returns
        -------
        named : dict
            A dictonary with the keys from `HEADERS` and values suitable for
            writing to the report.
        """

        protein_species = self.protein_species
        proteins = self.proteins
        if max_proteins is not None and len(protein_species) >= max_proteins:
            proteins = []
            protein_species = []

        return {
            'Group': self.group,
            'Rank': self.rank,
            'Type': self.type,
            'IFE id': self.ife_id,
            'PDB': self.pdb_id,
            'Chains': ', '.join(self.names),
            'Protein Species': ', '.join(str(p) for p in protein_species),
            'Protein Compound': ', '.join(str(p) for p in proteins),
            'RNA Species': self.species,
            'Nucleic Acid Compound': self.compound,
            'Name': '',
            'Observed Length': self.observed,
            'Experimental Length': self.experimental,
            'Experimental Sequence': self.sequence,
        }


class Info(object):
    def __init__(self, maker):
        self.session = maker

    @property
    def pdbs(self, **kwargs):
        if hasattr(self, '_pdb'):
            return self._pdb

        with self.session() as session:
            query = session.query(mod.ChainInfo).\
                order_by(mod.ChainInfo.pdb_id)

            results = it.imap(row2dict, query)
            grouped = it.groupby(results, op.itemgetter('pdb_id'))
            self._pdb = {}
            for pdb, group in grouped:
                chains = list(group)
                proteins = [c['chain_id'] for c in chains if is_protein(c)]
                rna = [c['chain_id'] for c in chains if is_rna(c)]
                self._pdb[pdb] = PdbInfo(proteins, rna)
        return self._pdb

    @property
    def chains(self, **kwargs):
        if hasattr(self, '_chains'):
            return self._chains

        self._chains = {}
        with self.session() as session:
            chain = mod.ChainInfo
            species = mod.ChainSpecies
            mapping = mod.SpeciesMapping
            exp_map = mod.ExpSeqChainMapping
            exp = mod.ExpSeqInfo
            units = mod.UnitInfo

            query = session.query(chain.chain_name,
                                  chain.chain_id,
                                  mapping.species_name.label('species'),
                                  chain.compound,
                                  chain.sequence,
                                  func.count(units.unit_id).label('observed'),
                                  exp.length.label('experimental'),
                                  ).\
                outerjoin(species, species.chain_id == chain.chain_id).\
                outerjoin(mapping, mapping.species_id == species.species_id).\
                outerjoin(exp_map, exp_map.chain_id == chain.chain_id).\
                outerjoin(exp, exp.exp_seq_id == exp_map.exp_seq_id).\
                outerjoin(units,
                          (units.pdb_id == chain.pdb_id) &
                          (units.chain == chain.chain_name) &
                          (units.unit_type_id.in_(['rna', 'aa']))).\
                group_by(chain.chain_id)

            for result in query:
                info = ChainInfo(result.chain_name,
                                 result.species,
                                 result.compound,
                                 result.observed,
                                 result.experimental,
                                 result.sequence)
                self._chains[result.chain_id] = info
        return self._chains


def sort_groups(data):
    key = op.attrgetter('group')
    grouped = it.groupby(sorted(data, key=key), key)
    ordering = []
    mapping = {}
    for name, members in grouped:
        members = sorted(members, key=op.attrgetter('rank'))
        mapping[name] = members
        rep_size = members[0].experimental
        ordering.append((name, rep_size))

    ordering.sort(key=op.itemgetter(1))
    ordered = []
    for name, _ in ordering:
        ordered.extend(mapping[name])
    return ordered


def groups(maker, release, resolution, **kwargs):
    """Create a report about the NR set.

    Parameters
    ----------
    maker : Session
        The session to use

    release : str
        The NR release to report on.

    resolution : str
        The resolution to use.

    Returns
    -------
    rows : list
        A list of dicts to write for the report.
    """

    info = Info(maker)
    with maker() as session:
        classes = mod.NrClasses
        chains = mod.NrChains
        ife = mod.IfeInfo
        ife_chains = mod.IfeChains
        query = session.query(classes.name,
                              chains.rank,
                              chains.ife_id,
                              ife.pdb_id,
                              ife_chains.chain_id,
                              ife_chains.index,
                              ).\
            join(chains, chains.nr_class_id == classes.nr_class_id).\
            join(ife, ife.ife_id == chains.ife_id).\
            join(ife_chains, ife.ife_id == ife_chains.ife_id).\
            filter(classes.nr_release_id == release).\
            filter(classes.resolution == resolution).\
            order_by(ife.ife_id)

        data = []
        for ife_id, group in it.groupby(query, op.attrgetter('ife_id')):
            chains = [row2dict(g) for g in group]
            result = chains[0]
            name = op.itemgetter('chain_id', 'index')
            chain_ids = set(name(c) for c in chains)
            chain_ids = sorted(chain_ids, key=op.itemgetter(1))
            chain_ids = [n[0] for n in chain_ids]
            current = Entry(result['name'],
                            result['rank'],
                            result['ife_id'],
                            result['pdb_id'],
                            chain_ids,
                            )
            current.add_info(info)
            data.append(current)
    return Report(GROUP_HEADERS, [d.named() for d in sort_groups(data)])


def pairs(maker, release, resolution, **kwargs):
    """Create a report about pairs in the NR set.
    """
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
    with maker() as session:
        query = session.query(classes.name.label('Group'),
                              classes.nr_release_id.label('Release'),
                              chains1.ife_id.label('IFE1'),
                              chains2.ife_id.label('IFE2'),
                              ccs.discrepancy.label('Discrepancy'),
                              corr.match_count.label('forward_match'),
                              rev_corr.match_count.label('rev_match'),
                              exp1.length.label('len1'),
                              exp2.length.label('len2'),
                              ).\
            join(chains1, chains1.nr_class_id == classes.nr_class_id).\
            join(chains2, chains2.nr_class_id == classes.nr_class_id).\
            join(ife1,
                 (ife1.ife_id == chains1.ife_id) &
                 (ife1.index == 0)
                 ).\
            join(ife2,
                 (ife2.ife_id == chains2.ife_id) &
                 (ife2.index == 0)
                 ).\
            join(mapping1, mapping1.chain_id == ife1.chain_id).\
            join(mapping2, mapping2.chain_id == ife2.chain_id).\
            join(exp1, exp1.exp_seq_id == mapping1.exp_seq_id).\
            join(exp2, exp2.exp_seq_id == mapping2.exp_seq_id).\
            outerjoin(ccs,
                      (ccs.chain_id_1 == ife1.chain_id) &
                      (ccs.chain_id_2 == ife2.chain_id),
                      ).\
            outerjoin(corr,
                      (corr.exp_seq_id_1 == mapping1.exp_seq_id) &
                      (corr.exp_seq_id_2 == mapping2.exp_seq_id),
                      ).\
            outerjoin(rev_corr,
                      (rev_corr.exp_seq_id_1 == mapping1.exp_seq_id) &
                      (rev_corr.exp_seq_id_2 == mapping2.exp_seq_id),
                      ).\
            filter(chains1.nr_chain_id != chains2.nr_chain_id).\
            filter(ife1.chain_id != ife2.chain_id).\
            filter(classes.resolution == resolution).\
            filter(classes.nr_release_id == release).\
            distinct().\
            order_by(classes.name, ife1.ife_id, ife2.ife_id)

        data = []
        for result in query:
            entry = row2dict(result)
            for_match = entry.pop('forward_match')
            rev_match = entry.pop('rev_match')
            len1 = entry.pop('len1')
            len2 = entry.pop('len2')
            entry['Alignment'] = None
            if for_match or rev_match:
                match = float(for_match or rev_match)
                entry['Alignment'] = match / float(min(len1, len2))
            data.append(entry)

    return Report(PAIRS_HEADERS, data)
