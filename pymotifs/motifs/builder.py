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
from pymotifs.utils.naming import ChangeCounter
from pymotifs.constants import MOTIF_GROUP_NAME

from pymotifs.motifs.cluster import ClusterMotifs


class BaseParser(core.Base):
    """A base parser for all motif csv files. This is callable and will produce
    a generator of dictonaries. If the row has a 'name' property we will yield
    the name and dictonary, otherwise we just yield the dictonary. This will
    also convert the field values specified by convert property.
    """

    """Name of the file to parser."""
    filename = None

    """List of fieldnames to use."""
    header = []

    """Converters, if any, to apply to fields in the file."""
    convert = {}

    def __convert__(self, entry):
        """Converts the values as requested.

        :param dict entry: The dictonary to convert.
        """

        data = {}
        for key, value in entry.items():
            converter = self.convert.get(key, str)
            data[key] = converter(value)
        return data

    def __call__(self, directory):
        """Parse the file in the given directory.

        :param str directory: The directory to get the file in.
        :yields: Each row in the file, after converting.
        """

        with open(os.path.join(directory, self.filename)) as raw:
            reader = csv.DictReader(raw, fieldnames=self.header)
            seen = False
            for row in reader:
                converted = self.__convert__(row)
                if 'name' in converted:
                    name = converted.pop('name')
                    yield name, converted
                else:
                    yield converted
                seen = True

            if not seen:
                raise core.InvalidState("No rows for %s" % self.filename)


class MotifListLoader(BaseParser):
    """Parse the MotifList.csv file.
    """
    filename = 'MotifList.csv'
    header = ['id', 'name']


class MotifPositionLoader(BaseParser):
    """Parse the MotifPositionLoader file.
    """
    filename = 'MotifPositions.csv'
    header = ['name', 'loop_id', 'unit_id', 'position']
    convert = {'position': int}


class MotifLoopOrder(BaseParser):
    """Parse the MotifLoopOrder file
    """
    filename = 'MotifLoopOrder.csv'
    header = ['name', 'loop_id', 'original_order', 'similarity_order']
    convert = {'original_order': int, 'similarity_order': int}


class MutualDiscrepancyLoader(BaseParser):
    """Parse the MutualDiscrepancy file.
    """
    filename = 'MutualDiscrepancy.csv'
    header = ['loop_id_1', 'discrepancy', 'loop_id_2']
    convert = {'discrepancy': float}


class BpSignaturesLoader(BaseParser):
    """Parse the MotifBpSignatures file.
    """
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
            '2d': None,
        }

    def loops(self, directory, data):
        """Load all loop to motif assignments in the given directory.
        """

        loader = MotifListLoader(self.config, self.session)
        c = 0
        name_dict = {}
        for name, loop in loader(directory):
            data[name]['members'].append(loop)
            c += 1
            name_dict[name] = 1
        self.logger.info('Found %d loops in MotifListLoader' % c)
        self.logger.info('Found %d names in MotifListLoader' % len(name_dict.keys()))

        return data

    def positions(self, directory, data):
        """Load all loop position information.
        """

        loader = MotifPositionLoader(self.config, self.session)
        c = 0
        name_dict = {}
        for name, entry in loader(directory):
            data[name]['positions'].append(entry)
            c += 1
            name_dict[name] = 1
        self.logger.info('Found %d positions in MotifPositionLoader' % c)
        self.logger.info('Found %d names in MotifPositionLoader' % len(name_dict.keys()))
        return data

    def ordering(self, directory, data):

        loader = MotifLoopOrder(self.config, self.session)
        c = 0
        name_dict = {}
        for name, entry in loader(directory):
            data[name]['ordering'].append(entry)
            c += 1
            name_dict[name] = 1
        self.logger.info('Found %d positions in MotifLoopOrder' % c)
        self.logger.info('Found %d names in MotifLoopOrder' % len(name_dict.keys()))
        return data

    def signature(self, directory, data):
        loader = BpSignaturesLoader(self.config, self.session)
        c = 0
        name_dict = {}
        for name, entry in loader(directory):
            data[name]['signature'] = entry['bp_signature']
            c += 1
            name_dict[name] = 1
        self.logger.info('Found %d signatures in BpSignaturesLoader' % c)
        self.logger.info('Found %d names in BpSignaturesLoader' % len(name_dict.keys()))
        return data

    def secondary_structures(self, directory, data):
        for name, motif in data.items():
            motif['2d'] = os.path.join(directory, '2ds', name + '.png')
        return data

    def __call__(self, release_id, directory):
        """Load and merge all motif information for the given release.
        """

        data = coll.defaultdict(ft.partial(self.empty_motif, release_id))
        data = self.loops(directory, data)
        data = self.positions(directory, data)
        data = self.ordering(directory, data)
        data = self.signature(directory, data)
        data = self.secondary_structures(directory, data)
        return data


class Known(core.Base):
    """A class to help load information from the database.
    """

    def loops(self, loop_type, ml_release_id):
        """Get a set of all known loops that were part of the given motif
        release.
        """
        motifs = self.motifs(loop_type, ml_release_id)
        loops = set()
        for motif in motifs:
            loops.update(motif.members)
        return loops

    def handles(self):
        """Get all handles that hav ever been used.
        """

        with self.session() as session:
            query = session.query(mod.MlMotifsInfo.handle).distinct()
            return set(result.handle for result in query)

    def names(self, loop_type, ml_release_id):
        """Get the names of motifs of the given type in the given release.
        """

        motifs = self.motifs(loop_type, ml_release_id)
        return set(m['name'] for m in motifs)

    def motifs(self, loop_type, ml_release_id):
        """Get a listing of all motifs in the given release.
        """

        with self.session() as session:
            motifs = mod.MlMotifsInfo
            query = session.query(
                motifs.handle,
                motifs.version,
                motifs.motif_id,
                mod.MlLoops.loop_id,
            ).\
                join(mod.MlLoops,
                     (mod.MlLoops.motif_id == motifs.motif_id) &
                     (mod.MlLoops.ml_release_id == motifs.ml_release_id)).\
                join(mod.LoopInfo,
                     mod.LoopInfo.loop_id == mod.MlLoops.loop_id).\
                filter(mod.MlLoops.ml_release_id == ml_release_id).\
                filter(mod.LoopInfo.type == loop_type)

            data = coll.defaultdict(lambda: {'name': None, 'members': []})
            for result in query:
                data[result.motif_id]['name'] = {
                    'full': result.motif_id,
                    'handle': result.handle,
                    'version': int(result.version),
                }
                data[result.motif_id]['members'].append({'id': result.loop_id})
            return sorted(data.values(), key=lambda d: d['name']['full'])


class Builder(core.Base):
    """This will build the data structures for a motif release. This will load
    and merge all motif data as well as compute names for the given
    """

    def motifs(self, loop_type, parent_id, release_id, directory):
        """Load all motif data and then name them.

        :param str parent_id: The id of the parent release.
        :param str release_id: The current release id.
        :param str directory: The directory all results are stored in.
        :returns: A list of the named motifs.
        """

        combiner = Combiner(self.config, self.session)
        known = Known(self.config, self.session)
        motifs = combiner(release_id, directory).values()
        parents = known.motifs(loop_type, parent_id)
        handles = known.handles()
        namer = Namer(self.config, self.session)
        named = namer(motifs, parents, handles)
        generate = MOTIF_GROUP_NAME.format
        for motif in named:
            motif['motif_id'] = generate(type=loop_type,
                                         handle=motif['name']['handle'],
                                         version=motif['name']['version'])
        return named

    def mutual_discrepancy(self, directory):
        """Load all mutual discrepancy information.

        :param str dictonary: The location of the clustering results.
        :returns: A list of all mutual discrepancies.
        """

        loader = MutualDiscrepancyLoader(self.config, self.session)
        return list(loader(directory))

    def parent_counts(self, parents, motifs):
        """Compute the number of changes between this release and the parent
        release.

        :param list parents: The parent motifs.
        :param list motifs: The motifs in the currrent release.
        :returns: A dictonary with entries for loops, motifs, and pdbs.
        """

        def as_pdb(group):
            return []

        counter = ChangeCounter(self.config, self.session)
        counts = counter(motifs, parents, pdbs=as_pdb)
        counts['loops'] = counts.pop('members')
        counts['motifs'] = counts.pop('groups')
        del counts['loops']['unchanged']
        del counts['pdbs']['unchanged']

        return counts

    def graph(self, directory):
        return os.path.join(directory, 'Supergroups.graphml')

    def __call__(self, loop_type, parent_id, release_id, loops,
                 directory=None):
        """Build the data structures for a motif release. This will load and
        name all motifs as well as load the mutual discrepancy information.

        :param str parent_id: The id of the parent release.
        :param str release_id: The current release id.
        :param str directory: The directory all results are stored in. If this
        is passed then the results there will be loaded, otherwise the motifs
        will be computed and placed in a directory.
        :returns: A dictonary with 2 keys, 'motifs' for a list of all motif
        from the release, and 'discrepancies' a list of all loop-loop
        discrepancies.
        """

        if not loops and directory is None:
            raise ValueError("Cannot build motifs without loops")

        if not directory:
            cluster = ClusterMotifs(self.config, self.session)
            directory = cluster(loop_type, loops, release_id)

        motifs = self.motifs(loop_type, parent_id, release_id, directory)
        known = Known(self.config, self.session)
        parents = known.motifs(loop_type, parent_id)

        return {
            'motifs': motifs,
            'loop_type': loop_type,
            'parent_counts': self.parent_counts(parents, motifs),
            'discrepancies': self.mutual_discrepancy(directory),
            'release': release_id,
            'parent': parent_id,
            'graph': self.graph(directory),
        }
