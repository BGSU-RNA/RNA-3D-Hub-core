"""Run all exporting functionality. This contains the:

export.cif_atom
    Write cifatom files.
export.interactions
    Write the interaction CSV
export.loops
    Write the loops CSV
"""

from pymotifs.core import StageContainer

from pymotifs.export.cifatom import Exporter as CifAtom
from pymotifs.export.interactions import Exporter as InteractionExporter
from pymotifs.export.loops import Exporter as LoopExporter
from pymotifs.export.pickle_units import Exporter as PickleExporter


class Exporter(StageContainer):
    """The actual stage to run."""

    stages = set([CifAtom, InteractionExporter, LoopExporter, PickleExporter])
