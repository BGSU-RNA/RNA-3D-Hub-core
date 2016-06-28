"""This module contains code to parse all the files that the motif pipeline
writes. It also contains the logic for naming new motifs using the parent
release.
"""

import os
import csv
import functools as ft
import collections as coll

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils.naming import Namer


class BaseParser(core.Base):
    filename = None
    header = []
    convert = {}

    def __convert__(self, entry):
        data = {}
        for key, value in entry.items():
            converter = self.convert.get(key, str)
            data[key] = converter(value)
        return data

    def __call__(self, directory):
        with open(os.path.join(directory, self.filename)) as raw:
            reader = csv.DictReader(raw, fieldnames=self.header)
            for row in reader:
                converted = self.__convert__(row)
                if 'name' in converted:
                    name = converted.pop('name')
                    yield name, converted
                else:
                    yield converted


class MotifListLoader(BaseParser):
    filename = 'MotifList.csv'
    header = ['id', 'name']


class MotifPositionLoader(BaseParser):
    filename = 'MotifPositions.csv'
    header = ['name', 'loop_id', 'unit_id', 'position']
    convert = {'position': int}


class MotifLoopOrder(BaseParser):
    filename = 'MotifLoopOrder.csv'
    header = ['name', 'loop_id', 'original_order', 'similarity_order']
    convert = {'original_order': int, 'similarity_order': int}


class MutualDiscrepancyLoader(BaseParser):
    filename = 'MutualDiscrepancy.csv'
    header = ['loop_id_1', 'discrepancy', 'loop_id_2']
    convert = {'discrepancy': float}


class BpSignaturesLoader(BaseParser):
    filename = 'MotifBpSignatures.csv'
    header = ['name', 'bp_signature']


class Combiner(core.Base):
    '''This is a class to combine all motif data into one single python data
    structure. This makes future work on doing things like naming and storing a
    lot easier.
    '''

    def empty_motif(self, release_id):
        return {
            'release_id': release_id,
            'members': [],
            'positions': [],
            'ordering': [],
            'signature': None,
        }

    def loops(self, directory, data):
        """Load all loop to motif assignments in the given directory.
        """

        loader = MotifListLoader(self.config, self.session)
        for name, loop in loader(directory):
            data[name]['members'].append(loop)
        return data

    def positions(self, directory, data):
        loader = MotifPositionLoader(self.config, self.session)
        for name, entry in loader(directory):
            data[name]['positions'].append(entry)
        return data

    def ordering(self, directory, data):
        loader = MotifLoopOrder(self.config, self.session)
        for name, entry in loader(directory):
            data[name]['ordering'].append(entry)
        return data

    def signature(self, directory, data):
        loader = BpSignaturesLoader(self.config, self.session)
        for name, entry in loader(directory):
            data[name]['signature'] = entry['bp_signature']
        return data

    def __call__(self, release_id, directory):
        data = coll.defaultdict(ft.partial(self.empty_motif, release_id))
        data = self.loops(directory, data)
        data = self.positions(directory, data)
        data = self.ordering(directory, data)
        data = self.signature(directory, data)
        return data


class Builder(core.Base):
    def load_motifs(self, release_id):
        with self.session() as session:
            query = session.query(mod.MlLoops).\
                filter_by(ml_release_id=release_id)

            data = coll.defaultdict(list)
            for result in query:
                data[result.motif_id].append({'id': result.loop_id})
            return data.values()

    def known_handles(self):
        with self.session() as session:
            query = session.query(mod.MlMotifInfo.handle).distinct()
            return set(result.handle for result in query)

    def naming(self, parents, groups, handles):
        namer = Namer(self.config, self.session)
        return namer(groups, parents, handles)

    def __call__(self, parent_id, release_id, directory):
        combiner = Combiner(self.config, self.session)
        motifs = combiner(release_id, directory)
        parents = self.load_motifs(parent_id)
        handles = self.known_handles()
        named = self.naming(parents, motifs, handles)
        return named
