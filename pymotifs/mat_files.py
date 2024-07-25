"""Generate requried mat files. This is a stage to generate the mat files that
are required for several other stages. We do this so we can force recomputing
of matlab files as well as make it easier to reason about what things require
matlab. All things that use matlab should depend on this stage.
"""
import os

from pymotifs import core
from pymotifs.utils import matlab
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
        filename = self.filename(pdb)
        if os.path.exists(filename):
            os.remove(filename)

    def data(self, pdb, **kwargs):

        # prevent trying to create Matlab interactions when filling in DNA releases
        nr_molecule_parent_current = kwargs.get('nr_molecule_parent_current','')
        if nr_molecule_parent_current and 'dna' in nr_molecule_parent_current.lower():
            raise core.Skip("no need to run Matlab on DNA structures")

        mlab = matlab.Matlab(self.config['locations']['fr3d_root'])
        mlab.zAddNTData(pdb)
        return None
