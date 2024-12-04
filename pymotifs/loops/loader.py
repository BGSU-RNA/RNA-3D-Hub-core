from pymotifs import core

# from pymotifs.loops.extractor import Loader as Extractor
# from pymotifs.loops.positions import Loader as PositionLoader
from pymotifs.loops.save_loops import Loader as Saver
from pymotifs.loops.quality import Loader as QALoader
# from pymotifs.loops.add_modified_nucleotides import Loader as ModifiedNucleotidesLoader


class Loader(core.StageContainer):
    # stages = set([Extractor, PositionLoader, QALoader])
    stages = set([Saver, QALoader])
