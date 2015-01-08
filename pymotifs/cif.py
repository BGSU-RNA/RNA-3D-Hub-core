"""
This module deals with cif files. This servers to provide a stage which loads
up cif files. This is useful because several other stages process cif files. By
declaring the stage her as their input the cif file need only be parsed once.
"""

import core
import utils

from download import Downloader

from fr3d.cif.reader import Cif as Reader


class CifStructure(core.Stage):
    depends_on = (Downloader)
    update_gap = False

    def __init__(self, *args):
        super(CifStructure, self).__init__(*args)
        self.filename = utils.CifFileFinder(self.config)

    def process(self, pdb):
        with open(self.filename(pdb), 'rb') as raw:
            return Reader(raw).structure()
