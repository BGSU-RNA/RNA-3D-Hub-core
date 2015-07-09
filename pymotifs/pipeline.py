from pymotifs import core

from pymotifs.downloader import Downloader
from pymotifs.export.cif_atom import Exporter as CifAtomExpoter
from pymotifs.pdbs import Loader as PdbLoader
from pymotifs.nts import Loader as NtLoader
from pymotifs.interactions import Loader as InterLoader
from pymotifs.chains import Loader as ChainLoader
from pymotifs.exp_seq import Loader as ExpSeqLoader
from pymotifs.correspondence import Loader as CorrLoader
from pymotifs.export import Exporter as GeneralExporter


class Loader(core.MultiLoader):
    dependencies = set([Downloader, CifAtomExpoter, PdbLoader, NtLoader,
                        InterLoader, ExpSeqLoader, ChainLoader, CorrLoader,
                        GeneralExporter])
