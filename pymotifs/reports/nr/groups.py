"""This module contains the logic to create the NR reports.
"""

import operator as op
import itertools as it
import collections as coll

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from sqlalchemy.sql.expression import func

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
    release : str
        The nr release id used.
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
    pdb_validation : dict
        A dict of structure level quality data.
    """

    def __init__(self, group, release, rank, ife_id, pdb_id, chain_ids):
        """Create a new `Entity`.

        Parameters
        ----------
        group : str
            The name of the equivlance class this
        release : str
            The NR release id used.
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
        self.release = release
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
        self.bp = 0.0
        self.nt = 0.0
        self.pdb_validation = {}

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
        self.pdb_validation = info.validation.get(self.pdb_id, {})

        ife_info = info.interactions[self.ife_id]
        self.bp = ife_info['bp']

        for index, chain_id in enumerate(self.chains):
            chain_info = info.chains[chain_id]
            self.names.append(chain_info.name)
            if index == 0:
                self.species = chain_info.species
                self.compound = chain_info.compound
                self.observed = chain_info.observed
                self.experimental = chain_info.experimental
                self.sequence = chain_info.sequence

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

        bp_nt = 0
        if self.observed:
            bp_nt = round(float(self.bp) / float(self.observed), 3)

        validation = lambda k: self.pdb_validation.get(k, None)

        return {
            'Group': self.group,
            'Release': self.release,
            'IFE id': self.ife_id,
            'BP/NT': bp_nt,
            'PDB': self.pdb_id,
            'Chains': ', '.join(self.names),
            'Protein Species': ', '.join(str(p) for p in protein_species),
            'Protein Compound': ', '.join(str(p) for p in proteins),
            'RNA Species': self.species,
            'Nucleic Acid Compound': self.compound,
            'Observed Length': self.observed,
            'Experimental Length': self.experimental,
            'Experimental Sequence': self.sequence,
            'Percent RSRZ Outliers': validation('percent_rsrz_outliers'),
            'Clashscore': validation('clashscore'),
            'Percent RNA Backbone Outliers': validation('percent_rota_outliers'),
            'Relative Percentile RSRZ': validation('relative_percentile_percent_rsrz_outliers'),
            'Relative Percentile Clashscore': validation('relative_percentile_clashscore'),
            'Relative Percentile Backbone Outliers': validation('relative_percentile_percent_rota_outliers'),
        }


class Info(object):
    def __init__(self, maker):
        self.session = maker

    @property
    def interactions(self):
        if hasattr(self, '_interactions'):
            return self._interactions

        with self.session() as session:
            query = session.query(mod.IfeInfo.ife_id,
                                  mod.IfeInfo.bp_count)

            self._interactions = {}
            for result in query:
                self._interactions[result.ife_id] = {
                    'bp': result.bp_count,
                }
        return self._interactions

    @property
    def pdbs(self):
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
    def validation(self):
        if hasattr(self, '_validation'):
            return self._validation

        self._validation = {}
        with self.session() as session:
            query = session.query(mod.PdbQuality)
            for result in query:
                entry = row2dict(result)
                pdb_id = entry.pop('pdb_id')
                self._validation[pdb_id] = entry
        return self._validation

    @property
    def chains(self):
        if not hasattr(self, '_chains'):
            self._chains = self.__chains__()
        return self._chains

    def __chains__(self):
        chains = {}
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
                chains[result.chain_id] = info
        return chains


class Groups(core.Reporter):
    headers = [
        'Original Index',
        'Group',
        'Release',
        'IFE id',
        'BP/NT',
        'PDB',
        'Chains',
        'Protein Species',
        'Protein Compound',
        'RNA Species',
        'Nucleic Acid Compound',
        'Observed Length',
        'Experimental Length',
        'Experimental Sequence',
        'Percent RSRZ Outliers',
        'Clashscore',
        'Percent RNA Backbone Outliers',
        'Relative Percentile RSRZ',
        'Relative Percentile Clashscore',
        'Relative Percentile Backbone Outliers'
    ]

    def sort_groups(self, data):
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

    def data(self, entry, **kwargs):
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

        release, resolution = entry
        info = Info(self.session)
        with self.session() as session:
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
                                release,
                                result['rank'],
                                result['ife_id'],
                                result['pdb_id'],
                                chain_ids,
                                )
                current.add_info(info)
                data.append(current)

        result = []
        for index, entry in enumerate(self.sort_groups(data)):
            entry = entry.named()
            entry['Original Index'] = index
            result.append(entry)
        return result
