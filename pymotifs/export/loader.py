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
from pymotifs.export.pickle_units_rna import Exporter as PickleURExporter
from pymotifs.export.pickle_pairs_rna import Exporter as PicklePRExporter
from pymotifs.export.pickle_unit_annotations import Exporter as PickleUAExporter
from pymotifs.export.ife_discrepancy import Exporter as IFEDiscrepancyExporter
from pymotifs.export.NA_datafile import Exporter as NA_datafileExporter

class Exporter(StageContainer):
    """The actual stage to run."""

#    stages = set([CifAtom, InteractionExporter, LoopExporter, PickleURExporter, PicklePRExporter])

    stages = set([CifAtom, InteractionExporter, LoopExporter, PickleURExporter, PicklePRExporter, PickleUAExporter, IFEDiscrepancyExporter, NA_datafileExporter])
