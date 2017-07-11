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


class Groups(core.Reporter):
    headers = [
        'Original Index',
        'Release',
        'Group',
        'PDB',
        'Title',
        'Chains',
        'Resolution',
        'Method',
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
            grouped = it.imap(row2dict, query)
            grouped = it.groupby(grouped, op.itemgetter('pdb_id'))
            for pdb_id, results in grouped:
                proteins = coll.defaultdict(list)
                for result in results:
                    for key, value in result.items():
                        if key != 'pdb_id':
                            proteins[key].append(value)

                proteins = dict(proteins)
                for key, value in proteins.items():
                    proteins[key] = ', '.join(value)
                data[pdb_id] = proteins
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
                mod.PdbInfo.resolution.label('Resolution'),
                mod.PdbInfo.experimental_technique.label('Method'),
                mod.PdbInfo.title.label('Title'),
            ).filter(mod.PdbInfo.pdb_id.in_(pdb_ids))

            data = {}
            for result in query:
                entry = row2dict(result)
                method = entry.pop('Method')
                if method == 'X-RAY DIFFRACTION':
                    method = 'x-ray'
                elif method == 'SOLUTION NMR':
                    method = 'nmr'
                elif method == 'ELECTRON MICROSCOPY':
                    method = 'cyro-em'
                else:
                    method = method
                entry['Method'] = method
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
