"""This is a loader that will copy the annotation data from a parent into the
child as needed
"""

from datetime import datetime

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.motifs.utils import BaseLoader
from pymotifs.motifs.info import Loader as InfoLoader
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    """
    Loader for the motif annotations.
    """

    dependencies = set([ReleaseLoader, InfoLoader])

    @property
    def table(self):
        return mod.MlMotifAnnotations

    def has_data(self, *args, **kwargs):
        return False

    def remove(self, *args, **kwargs):
        logger.exception("Cannot easily clean up motif annotations")

    def known(self):
        annotations = {}
        with self.session() as session:
            query = session.query(self.table).\
                order_by(self.table.date.asc())
            for result in query:
                current = row2dict(result)
                for key, value in current.items():
                    if value == '':
                        current[key] = None
                annotations[result.motif_id] = current
        return annotations

    def parent_annotations(self, known, motif, name):
        values = []
        for parent in motif['parents']:
            parent_id = parent['name']['full']
            if parent_id not in known:
                msg = "Parent (%s) of %s missing annotations"
                raise core.InvalidState(msg % (parent_id, motif_id))
            current = known[parent_id]
            value = current[name]
            if value:
                values.extend(value.split('|'))

        values = [v.strip() for v in values if v]
        if len(values) != len(motif['parents']):
            self.logger.info(
                "Not all parents have '%s' annotations for %s, not propagating",
                name,
                motif['motif_id'],
            )
            return None

        values = set(values)
        if len(values) > 1:
            self.logger.info(
                "Inconsistent parent '%s' annotations for %s (%s), not propagating",
                name,
                motif['motif_id'],
                '|'.join(sorted(values)),
            )
            return None

        if not values:
            self.logger.info(
                "No parent %s annotation for %s",
                name,
                motif['motif_id'],
            )
            return None
        return values.pop()

    def new_annotations(self, known, cached):
        annotations = []
        for motif in cached['motifs']:
            motif_id = motif['motif_id']
            if motif_id in known:
                continue
            annotations.append({
                'motif_id': motif_id,
                'common_name': self.parent_annotations(known, motif, 'common_name'),
                'annotation': self.parent_annotations(known, motif, 'annotation'),
                'bp_signature': motif['signature'],
                'date': datetime.now(),
            })
        return annotations

    def data(self, pair, **kwargs):
        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No cached data")

        if cached['release'] == cached['parent']:
            raise core.Skip("No annotations for first release")

        known = self.known()
        if not known:
            raise core.Skip("No existing annotations")

        return self.new_annotations(known, cached)
