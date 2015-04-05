"""This is a class to generate the mat files that are required for several
other stages. We do this so we can force recomputing of matlab files as well as
make it easier to reason about what things require matlab. All things that use
matlab should depend on this stage.
"""
import os

from pymotifs import core
from pymotifs.export.cifatom import Exporter as CifAtom


class Loader(core.Loader):
    allow_no_data = True
    dependencies = set([CifAtom])

    def filename(self, pdb):
        return os.path.join(self.config['locations']['fr3d_root'],
                            'PrecomputedData', pdb + '.mat')

    def has_data(self, pdb, **kwargs):
        return os.path.exists(self.filename(pdb))

    def remove(self, pdb, **kwargs):
        os.remove(self.filename(pdb))

    def data(self, pdb, **kwargs):
        matlab = core.Matlab()
        matlab.zAddNTData(pdb)
        return None
