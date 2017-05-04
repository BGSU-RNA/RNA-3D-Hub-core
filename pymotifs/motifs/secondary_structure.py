"""
This is a loader which will move the generated motif files into the correct
directories for the website.
"""
import os
import shutil

from pymotifs import core
from pymotifs.motifs.release import Loader as ReleaseLoader


class Loader(core.Loader):
    allow_no_data = True
    dependencies = set([ReleaseLoader])

    def to_process(self, pdbs, **kwargs):
        """Compute the data to process. The input PDB's are ignored and instead
        the cache is examine for motif data, IL and HL to import. This will
        produce a list of tuples of the motifs to import.

        Parameters
        ----------
        pdbs : list
            Ingored

        Returns
        -------
        A list of tuples like [('IL', '1.0'), ('HL', '1.0'] to import.
        """
        current, _ = ReleaseLoader(self.config, self.session).current_id()
        data = []
        for loop_type in ReleaseLoader.types:
            cached = self.cached(loop_type)
            if not cached:
                raise core.InvalidState("No cached data")

            if cached['release'] != current:
                raise core.InvalidState("Caching does not match excepted ID")
            data.append((loop_type, current))
        return data

    def final_directory(self, loop_type, release):
        return os.path.join(
            self.config['locations']['2ds_destination'],
            loop_type + release
        )

    def final_location(self, release, motif):
        return os.path.join(
            self.final_directory(release, motif['type']),
            motif['motif_id'] + '.png',
        )

    def has_data(self, pair, **kwargs):
        return os.path.exists(self.final_directory(*pair))

    def remove(self, pair, **kwargs):
        filename = self.final_directory(*pair)
        if os.path.exists(filename):
            os.rmdir(filename)

    def data(self, pair, **kwargs):
        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No cached data")

        directory = self.final_directory(loop_type, release)
        if not directory:
            os.makedirs(directory)

        for motif in cached['motifs']:
            shutil.copy(motif['2d'], self.final_location(motif))
        return None