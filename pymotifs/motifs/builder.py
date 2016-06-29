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
        }

    def loops(self, directory, data):
        """Load all loop to motif assignments in the given directory.
        """

        loader = MotifListLoader(self.config, self.session)
        for name, loop in loader(directory):
            data[name]['members'].append(loop)
        return data

    def positions(self, directory, data):
        """Load all loop position information.
        """

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
        """Load and merge all motif information for the given release.
        """

        data = coll.defaultdict(ft.partial(self.empty_motif, release_id))
        data = self.loops(directory, data)
        data = self.positions(directory, data)
        data = self.ordering(directory, data)
        data = self.signature(directory, data)
        return data


class Builder(core.Base):
    """This will build the data structures for a motif release. This will load
    and merge all motif data as well as compute names for the given
    """

    def load_motifs(self, release_id):
        """Load all motifs for the given release id. This is used later to
        generate naming for the motifs.
        """

        with self.session() as session:
            query = session.query(mod.MlLoops).\
                filter_by(ml_release_id=release_id)

            data = coll.defaultdict(lambda: {'name': None, 'members': []})
            for result in query:
                data[result.motif_id]['name'] = result.motif_id
                data[result.motif_id]['members'].append({'id': result.loop_id})
            return data.values()

    def known_handles(self):
        """Get all known handles already in use. This is used in producing new
        names for the motifs.

        :returns: A set of all known handles.
        """

        with self.session() as session:
            query = session.query(mod.MlMotifInfo.handle).distinct()
            return set(result.handle for result in query)

    def motifs(self, parent_id, release_id, directory):
        """Load all motif data and then name them.

        :param str parent_id: The id of the parent release.
        :param str release_id: The current release id.
        :param str directory: The directory all results are stored in.
        :returns: A list of the named motifs.
        """

        combiner = Combiner(self.config, self.session)
        motifs = combiner(release_id, directory).values()
        parents = self.load_motifs(parent_id)
        handles = self.known_handles()
        namer = Namer(self.config, self.session)
        return namer(motifs, parents, handles)

    def mutual_discrepancy(self, directory):
        """Load all mutual discrepancy information.

        :param str dictonary: The location of the clustering results.
        :returns: A list of all mutual discrepancies.
        """

        loader = MutualDiscrepancyLoader(self.config, self.session)
        return list(loader(directory))

    def __call__(self, parent_id, release_id, directory):
        """Build the data structures for a motif release. This will load and
        name all motifs as well as load the mutual discrepancy information.

        :param str parent_id: The id of the parent release.
        :param str release_id: The current release id.
        :param str directory: The directory all results are stored in.
        :returns: A dictonary with 2 keys, 'motifs' for a list of all motif
        from the release, and 'discrepancies' a list of all loop-loop
        discrepancies.
        """

        return {
            'motifs': self.motifs(parent_id, release_id, directory),
            'discrepancies': self.mutual_discrepancy(directory),
        }
