"""Determine which interactions correspond in structures.
"""

from pymotifs import core
from pymotifs import models as mod

from pymotifs.models import UnitPairsInteractions as Ints
from pymotifs.models import CorrespondenceInteractions as CorrInts

from pymotifs.utils.correspondence import Helper


class Loader(core.Loader):

    def interactions(self, pdb):
        with self.session() as session:
            query = session.query(Ints.unit_pairs_interactions_id, Ints.unit_id_1, Ints.unit_id_2).\
                filter(Ints.pdb_id == pdb)

            data = []
            for result in query:
                data.append({
                    'id': result.id,
                    'unit1': result.unit_id_1,
                    'unit2': result.unit_id_2
                })
            return data

    def has_data(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(CorrInts).\
                join(Ints, Ints.id == CorrInts.interaction_id_1).\
                filter(Ints.pdb_id == pdb)

            return bool(query.count())

    def remove(self, pdb, **kwargs):
        with self.session() as session:
            query = session.query(Ints).\
                filter(Ints.pdb_id == pdb)
            ids = [result.id for result in query]

        with self.session() as session:
            session.query(CorrInts).\
                filter(CorrInts.correspondence_interactions_id.in_(ids)).\
                delete(synchronize_session=False)

    def compare(self, corr_id, first, second, mapping):
        known = {}
        for interaction in second:
            unit1 = interaction['unit1']
            unit2 = interaction['unit2']
            if unit1 in mapping and unit2 in mapping:
                known[(mapping[unit1], mapping[unit2])] = interaction

        for interaction in first:
            data = {
                'correspondence_id': corr_id,
                'interaction_id_1': interaction['id'],
                'interaction_id_2': None
            }

            key = (interaction['unit1'], interaction['unit2'])
            if key in known:
                data['interaction_id_2'] = known[key]['id']

            yield data

    def data(self, pdb1, **kwargs):

        util = Helper(self.session.maker)
        for pdb2 in util.pdbs(pdb1):
            corr_id, mapping = util.pdb_mapping(pdb1, pdb2)
            if not mapping:
                raise core.InvalidState("No mapping for %s %s", (pdb1, pdb2))

            interactions1 = self.interactions(pdb1)
            interactions2 = self.interactions(pdb2)
            for pair in self.compare(corr_id, interactions1, interactions2,
                                     mapping):
                yield CorrInts(**pair)
