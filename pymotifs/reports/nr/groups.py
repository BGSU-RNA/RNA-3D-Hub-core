"""This module contains the logic to create the NR reports.
"""

import operator as op
import itertools as it
import collections as coll
from copy import deepcopy

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.nr.representatives.using_quality import CompScore

import numpy as np

from sqlalchemy.sql.expression import func

compound = op.itemgetter('compound')
species = op.itemgetter('source')
macro_type = op.itemgetter('entity_macromolecule_type')
is_type = lambda t: lambda c: macro_type(c) == t
is_protein = is_type('Polypeptide(L)')
is_rna = is_type('Polyribonucleotide (RNA)')



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
        self.resolution = None

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
        self.resolution = pdb_info.resolution

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

    def validation(self, key):
        return self.pdb_validation.get(key, None)

    @property
    def validation_relative_mean(self):
        valid = self.validation
        measures = [
            'relative_percentile_percent_rsrz_outliers',
            'relative_percentile_clashscore',
            'relative_percentile_percent_rota_outliers',
        ]
        measures = [valid(m) for m in measures if valid(m) is not None]
        if measures:
            return np.mean(measures)
        return None

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

        relative_mean = self.validation_relative_mean
        if relative_mean is not None:
            relative_mean = round(relative_mean, 2)

        return {
            'Group': self.group,
            'Release': self.release,
            'IFE id': self.ife_id,
            'BP/NT': bp_nt,
            'PDB': self.pdb_id,
            'Resolution': self.resolution,
            'Chains': ', '.join(self.names),
            'Protein Species': ', '.join(str(p) for p in protein_species),
            'Protein Compound': ', '.join(str(p) for p in proteins),
            'RNA Species': self.species,
            'Nucleic Acid Compound': self.compound,
            'Observed Length': self.observed,
            'Experimental Length': self.experimental,
            'Experimental Sequence': self.sequence,
            'Percent RSRZ Outliers': self.validation('percent_rsrz_outliers'),
            'Clashscore': self.validation('clashscore'),
            'Percent RNA Backbone Outliers': self.validation('percent_rota_outliers'),
            'Relative Percentile RSRZ': self.validation('relative_percentile_percent_rsrz_outliers'),
            'Relative Percentile Clashscore': self.validation('relative_percentile_clashscore'),
            'Relative Percentile Backbone Outliers': self.validation('relative_percentile_percent_rota_outliers'),
            'Relative Percentile Mean': relative_mean,
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
            query = session.query(mod.ChainInfo.chain_id,
                                  mod.ChainInfo.entity_macromolecule_type,
                                  mod.ChainInfo.pdb_id,
                                  mod.PdbInfo.resolution).\
                join(mod.PdbInfo, mod.PdbInfo.pdb_id == mod.ChainInfo.pdb_id).\
                order_by(mod.ChainInfo.pdb_id)

            results = it.imap(row2dict, query)
            grouped = it.groupby(results, op.itemgetter('pdb_id'))
            self._pdb = {}
            for pdb, group in grouped:
                chains = list(group)
                proteins = [c['chain_id'] for c in chains if is_protein(c)]
                rna = [c['chain_id'] for c in chains if is_rna(c)]
                resolution = chains[0]['resolution']
                self._pdb[pdb] = PdbInfo(proteins, rna, resolution)
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
                                  func.count(units.unit_id.distinct()).label('observed'),
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
        'Release',
        'Group',
        'PDB',
        'Chains',
        'Resolution',
        'Observed Length',
        'Experimental Length',
        'Experimental Sequence',
        'BP/NT',
        'Nucleic Acid Compound',
        'RNA Species',
        'Protein Species',
        'Protein Compound',
        'Clashscore',
        'Compscore',
        'Average RSR',
        'Percent Clash',
        'Average RSCC',
        'Rfree',
    ]

    def class_property(self, ifes, name):
        return {ife[name] for ife in ifes}

    def quality_data(self, ifes):
        compscore = self._create(CompScore)
        members = deepcopy(ifes)
        compscore.load_quality(members)
        data = {}
        for ife in ifes:
            quality = ife['quality']
            data[ife['id']] = {
                'Clashscore': quality['clashscore'],
                'Compscore': compscore.compscore(ife),
                'Average RSR': quality['average_rsr'],
                'Percent Clash': quality['percent_clash'],
                'Average RSCC': quality['average_rscc'],
                'Rfree': quality['rfree'],
            }
        return data

    def ife_info(self, ifes):
        ife_ids = self.class_property(ifes, 'id')
        with self.session() as session:
            query = session.query(
                mod.IfeInfo.ife_id,
                mod.IfeInfo.bps,
                mod.IfeInfo.length.label('Observed Length')
            ).filter(mod.IfeInfo.ife_id.in_(ife_ids))

            data = {}
            for result in query:
                entry = row2dict(result)
                nt = entry['Observed Length']
                bp = entry.pop('bps')
                ife_id = entry.pop('ife_id')
                entry['BP/NT'] = float(bp) / float(nt)
                chain_ids = ife_id.split('+')
                chains = [p.split('|')[-1] for p in chain_ids]
                entry['Chains'] = ', '.join(chains)
                data[ife_id] = entry
        return data

    def protein_info(self, ifes):
        pdb_ids = self.class_property(ifes, 'pdb_id')
        with self.session() as session:
            query = session.query(
                mod.ChainInfo.pdb_id,
                mod.ChainInfo.compound.label('Protein Compound'),
                mod.SpeciesInfo.species_name.label('Protein Species'),
            ).\
                join(mod.ChainSpecies,
                     mod.ChainSpecies.chain_id == mod.ChainInfo.chain_id).\
                join(mod.SpeciesInfo,
                     mod.SpeciesInfo.species_id == mod.ChainSpecies.species_id).\
                filter(mod.ChainInfo.pdb_id.in_(pdb_ids)).\
                filter(mod.ChainInfo.entity_macromolecule_type == '')

            data = {}
            for result in query:
                entry = row2dict(result)
                pdb_id = entry.pop('pdb_id')
                data[pdb_id] = entry
            return data

    def chain_info(self, ifes):
        chain_ids = self.class_property(ifes, 'chain_id')
        with self.session() as session:
            query = session.query(
                mod.ChainInfo.chain_id,
                mod.ChainInfo.sequence.label('Experimental Sequence'),
                mod.ChainInfo.compound.label('Nucleic Acid Compound'),
                mod.SpeciesInfo.species_name.label('RNA Species'),
            ).\
                join(mod.ChainSpecies,
                     mod.ChainSpecies.chain_id == mod.ChainInfo.chain_id).\
                join(mod.SpeciesInfo,
                     mod.SpeciesInfo.species_id == mod.ChainSpecies.species_id).\
                filter(mod.ChainInfo.chain_id.in_(chain_ids))

            data = {}
            for result in query:
                entry = row2dict(result)
                entry['Experimental Length'] = len(entry['Experimental Sequence'])
                chain_id = entry.pop('chain_id')
                data[chain_id] = entry
        return data

    def pdb_info(self, ifes):
        pdb_ids = self.class_property(ifes, 'pdb_id')
        with self.session() as session:
            query = session.query(
                mod.PdbInfo.pdb_id.label('PDB'),
                mod.PdbInfo.resolution.label('Resolution')
            ).filter(mod.PdbInfo.pdb_id.in_(pdb_ids))

            data = {}
            for result in query:
                entry = row2dict(result)
                data[entry['PDB']] = entry
        return data

    def load_nr_classes(self, release, resolution):
        with self.session() as session:
            query = session.query(
                mod.NrChains.index,
                mod.NrChains.ife_id.label('id'),
                mod.IfeInfo.pdb_id,
                mod.IfeInfo.length,
                mod.IfeChains.chain_id,
                mod.NrClasses.name,
            ).\
                join(mod.IfeInfo, mod.IfeInfo.ife_id == mod.NrChains.ife_id).\
                join(mod.IfeChains,
                     mod.IfeChains.ife_id == mod.IfeInfo.ife_id).\
                join(mod.NrClasses,
                     mod.NrClasses.nr_class_id == mod.NrChains.nr_class_id).\
                filter(mod.NrClasses.nr_release_id == release).\
                filter(mod.NrClasses.resolution == resolution).\
                filter(mod.IfeChains.index == 0)

            data = []
            for result in query:
                entry = row2dict(result)
                entry['rep'] = (entry['index'] == 0)
                data.append(entry)
        return data

    def order_nr_classes(self, nr_classes):
        def key(nr_class):
            rep = next(ife for ife in nr_class if ife['rep'] is True)
            return rep['length']
        return sorted(nr_classes, key=key)

    def nr_classes(self, release, resolution):
        nr_classes = self.load_nr_classes(release, resolution)
        return self.order_nr_classes(nr_classes)

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
        nr_classes = self.nr_classes(release, resolution)
        index = 0
        for nr_class in nr_classes:
            pdb_info = self.pdb_info(nr_class)
            chain_info = self.chain_info(nr_class)
            ife_info = self.ife_info(nr_class)
            quality_data = self.quality_data(nr_class)
            protein_data = self.protein_data(nr_class)
            for ife in nr_class:
                data = {
                    'Original Index': index,
                    'Release': release,
                    'Group': ife['name'],
                }
                data.update(pdb_info[ife['pdb_id']])
                data.update(chain_info[ife['chain_id']])
                data.update(ife_info[ife['id']])
                data.update(quality_data[ife['id']])
                data.update(protein_data[ife['pdb_id']])
                yield data
                index + 1
