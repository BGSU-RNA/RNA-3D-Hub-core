import logging

import core
import utils
from models import ObsSequenceMapping as Mapping

from rnastructure.tertiary.cif import CIF

logger = logging.getLogger(__name__)


class Loader(core.Loader):
    name = 'experimental_mapping'
    update_gap = False
    insert_max = 5000

    def __init__(self, config, maker):
        self.finder = utils.CifFileFinder(config)
        super(Loader, self).__init__(config, maker)

    def has_data(self, pdb):
        with self.session() as session:
            query = session.query(Mapping).filter_by(pdb=pdb)
            return bool(query.count())

    def remove(self, pdb):
        with self.session() as session:
            session.query(Mapping).filter_by(pdb=pdb).delete()

    def data(self, pdb, **kwargs):
        mapping = []
        seen = set()
        with open(self.finder(pdb), 'rb') as raw:
            cif = CIF(raw)
            for chain in cif.chains():
                seq_mapping = chain.experimental_sequence_mapping()
                for (_, seq_id, unit_id) in seq_mapping:
                    if unit_id not in seen:
                        obj = Mapping(pdb=pdb, unit_id=unit_id,
                                      sequence_unit_id=seq_id)
                        mapping.append(obj)
                        seen.add(unit_id)
        return mapping
