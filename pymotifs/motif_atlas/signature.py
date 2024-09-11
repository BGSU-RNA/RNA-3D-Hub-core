"""This is a loader that will add basepair signatures as needed
"""

from datetime import datetime

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

from pymotifs.motif_atlas.utils import BaseLoader
from pymotifs.motif_atlas.info import Loader as InfoLoader
from pymotifs.motif_atlas.release import Loader as ReleaseLoader


class Loader(BaseLoader):
    """
    Loader for the motif annotations.
    """

    dependencies = set([ReleaseLoader, InfoLoader])

    """We allow for no data to be written when appropriate,
    for example, when all groups are the same as before"""
    allow_no_data = True

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


    def new_annotations(self, known, cached):
        annotations = []
        for motif in cached['motifs']:
            motif_id = motif['motif_id']
            if motif_id in known:
                continue
            annotations.append({
                'motif_id': motif_id,
                'common_name': '',
                'annotation': '',
                'bp_signature': motif['signature'],
                'date': datetime.now(),
            })
        return annotations

    def data(self, pair, **kwargs):
        loop_type, release = pair

        # retrieve data loaded from csv files
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No cached data")

        # query database for groups with signatures already
        known = self.known()

        return self.new_annotations(known, cached)
