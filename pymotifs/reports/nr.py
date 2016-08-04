import operator as op
import functools as ft
import itertools as it
import collections as coll

from pymotifs import models as mod
from pymotifs.utils import row2dict

from sqlalchemy.sql.expression import func


HEADERS = [
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

compound = op.itemgetter('compound')
species = op.itemgetter('source')
macro_type = op.itemgetter('macromolecule_type')
is_type = lambda t: lambda c: macro_type(c) == t
is_protein = is_type('Polypeptide(L)')
is_rna = is_type('Polyribonucleotide (RNA)')

PdbInfo = coll.namedtuple('PdbInfo', ['protein_chains', 'rna_chains'])
ChainInfo = coll.namedtuple('ChainInfo', ['name', 'species', 'compound',
                                          'observed', 'experimental',
                                          'sequence'])


class Entry(object):
    def __init__(self, rank, ife_id, pdb_id, chain_ids):
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

    def named(self, max_proteins=2):
        protein_species = self.protein_species
        proteins = self.proteins
        if len(protein_species) > max_proteins:
            proteins = []
            protein_species = []
        return {
            'Rank': self.rank,
            'Type': self.type,
            'IFE id': self.ife_id,
            'PDB': self.pdb_id,
            'Chains': ', '.join(self.names),
            'Protein Species': ', '.join(protein_species),
            'Protein Compound': ', '.join(proteins),
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
                                  # exp.normalized_length.label('observed'),
                                  exp.length.label('experimental'),
                                  ).\
                outerjoin(species, species.chain_id == chain.chain_id).\
                outerjoin(mapping, mapping.species_id == species.species_id).\
                outerjoin(exp_map, exp_map.chain_id == chain.chain_id).\
                outerjoin(exp, exp.exp_seq_id == exp_map.exp_seq_id).\
                join(units,
                     (units.pdb_id == chain.pdb_id) &
                     (units.chain == chain.chain_name)).\
                filter(units.unit_type_id == 'rna').\
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


def sorter(entry):
    pass


def report(maker, release, resolution, **kwargs):
    """Create a report about the NR set
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
            join(chains, chains.nr_chain_id == classes.nr_class_id).\
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
            current = Entry(result['rank'],
                            result['ife_id'],
                            result['pdb_id'],
                            chain_ids,
                            )
            current.add_info(info)
            data.append(current)
    return [d.named() for d in data]
