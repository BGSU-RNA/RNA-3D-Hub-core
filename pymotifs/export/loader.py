"""
Run all exporting functionality, including:

export.cif_atom
    Write cifatom files.
export.interactions
    Write the interaction CSV
export.loops
    Write the loops CSV
"""

from pymotifs.core import StageContainer

# from pymotifs.export.cifatom import Exporter as CifAtom
# from pymotifs.export.interactions import Exporter as InteractionExporter  # no longer needed by NDB, skip
# from pymotifs.export.loops import Exporter as LoopExporter   # would need to be updated
# from pymotifs.export.pickle_units_rna import Exporter as PickleURExporter
from pymotifs.export.pickle_pairs_na import Exporter as na_pairs_exporter
from pymotifs.export.pickle_unit_annotations import Exporter as PickleUAExporter
from pymotifs.export.ife_discrepancy import Exporter as IFEDiscrepancyExporter
from pymotifs.export.NA_datafile import Exporter as NA_datafileExporter
from pymotifs.export.pickle_units_na import Exporter as na_units_exporter
from pymotifs.export.pickle_units_protein import Exporter as protein_units_exporter

class Exporter(StageContainer):
    stages = set([na_pairs_exporter, PickleUAExporter, IFEDiscrepancyExporter, NA_datafileExporter, na_units_exporter, protein_units_exporter])
