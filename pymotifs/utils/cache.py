"""Some utilities for dealing with caching data. This focuses on caching cif
data and structure data.
"""

import os
import abc
import cPickle as pickle

from fr3d.cif.reader import Cif
from fr3d.cif import persist


class CacheData(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, config):
        self.config = config

    @abc.abstractmethod
    def cachename(self, pdb):
        pass

    @abc.abstractmethod
    def filename(self, pdb):
        pass

    @abc.abstractmethod
    def data(self, pdb):
        pass

    def clear(self):
        for filename in os.listdir(self.location):
            parts = os.path.splitext(filename)
            if parts[1] == '.pickle':
                os.unlink(os.path.join(self.location, filename))

    def cache(self, cache_file, data):
        with open(cache_file, 'w') as out:
            pickle.dump(data, out)

    def __call__(self, pdb):
        cache_file = self.cachename(pdb)
        data_file = self.filename(pdb)

        if not os.path.exists(cache_file) or \
                os.path.getmtime(cache_file) < os.path.getmtime(data_file):

            data = self.data(cache_file, pdb)
            self.cache(data)
            return data

        with open(cache_file, 'rb') as raw:
            return pickle.load(raw)


class CifData(CacheData):
    def __init__(self, config):
        self.cache = PickleFileFinder(config, strict=False)
        self.cif = CifFileFinder(config)

    def filename(self, pdb):
        return self.cif(pdb)

    def cache(self, filename, cif):
        with open(filename, 'rb') as out:
            persist.serialize(out, cif)

    def data(self, pdb):
        with open(self.name(pdb), 'rb') as raw:
            return Cif(raw)


class StructureData(CifData):
    def data(self, pdb):
        cif = super(CifData, self).data(pdb)
        structure = cif.structure()
        structure.infer_hydrogens()
        return structure
