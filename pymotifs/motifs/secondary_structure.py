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
        type = motif['motif_id'].split('_')[0]
        return os.path.join(
            self.final_directory(type, release),
            motif['motif_id'] + '.png',
        )

    def web_directory(self, loop_type, release):
        return os.path.join(
            self.config['locations']['web_dir'],
            loop_type + release
        )

    def web_location(self, release, motif):
        type = motif['motif_id'].split('_')[0]
        return os.path.join(
            self.web_directory(type, release),
            motif['motif_id'] + '.png',
        )

    def has_data(self, pair, **kwargs):
        return False

    def remove(self, pair, **kwargs):
        filename = self.final_directory(*pair)
        if os.path.exists(filename):
            os.rmdir(filename)

    def mark_processed(self, *args, **kwargs):
        return None

    def data(self, pair, **kwargs):
        loop_type, release = pair
        cached = self.cached(loop_type)
        if not cached:
            raise core.InvalidState("No cached data")

        directory = self.final_directory(loop_type, release)
        self.logger.debug("Use directory %s...", directory)
        if not os.path.isdir(directory):
            os.makedirs(directory)
            self.logger.debug("Needed to create directory %s", directory)

        if not os.path.isdir(directory):
            self.logger.debug("Failed to create directory %s", directory)

        web_dir = self.web_directory(loop_type, release)
        self.logger.debug("Web directory: %s", web_dir)
        if not os.path.isdir(web_dir):
            os.makedirs(web_dir)
            self.logger.debug("Needed to create directory %s", web_dir)

        if not os.path.isdir(web_dir):
            self.logger.debug("Failed to create directory %s", web_dir)

        for motif in cached['motifs']:
            location = self.final_location(release, motif)
            web_location = self.web_location(release, motif)
            self.logger.debug("Copying motif 2d %s %s to %s",
                              motif['motif_id'], motif['2d'], location)
            shutil.copy(motif['2d'], location)
            self.logger.debug("Copying motif 2d %s %s to %s",
                              motif['motif_id'], motif['2d'], web_location)
            shutil.copy(motif['2d'], web_location)

        return None
