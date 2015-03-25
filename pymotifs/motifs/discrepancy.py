from pymotifs.motifs.base import CsvLoader
from pymotifs.models import MlMutualDiscrepancy


class Loader(CsvLoader):
    name = None
    table = MlMutualDiscrepancy
    headers = ['loop_id1', 'discrepancy', 'loop_id2']
