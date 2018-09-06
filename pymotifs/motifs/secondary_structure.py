"""
This is a loader which will copy the generated motif files into the correct
directories for the website and JAR3D work.
"""
import os
import re
import shutil

from pprint import pformat

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

    def mat_location(self, mat_dir, release, motif):
        type = motif['motif_id'].split('_')[0]
        return os.path.join(
            mat_dir,
            motif['motif_id'] + '.mat',
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

        self.logger.debug("Cached data for %s: %s", loop_type, pformat(cached['motifs']))

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
            self.logger.debug("Needed to create web directory %s", web_dir)

        if not os.path.isdir(web_dir):
            self.logger.debug("Failed to create web directory %s", web_dir)

        mat_check = 0

        for motif in cached['motifs']:
            if mat_check == 0:
                mat_dir = re.sub('/2ds/.*$','/mat',motif['2d'])
                self.logger.debug("mat directory: %s", mat_dir)
                if not os.path.isdir(mat_dir):
                    os.makedirs(mat_dir)
                    self.logger.debug("Needed to create mat directory %s", mat_dir)
                    mat_check = 1

                if not os.path.isdir(mat_dir):
                    self.logger.debug("Failed to create mat directory %s", mat_dir)

            location = self.final_location(release, motif)
            web_location = self.web_location(release, motif)
            mat_file = motif['2d'].replace("2ds","Groups").replace("png","mat")
            mat_out = self.mat_location(mat_dir, release, motif)
            #mat_out = os.path.join(mat_dir, mat_file)
            self.logger.debug("Copying motif 2d %s %s to %s",
                              motif['motif_id'], motif['2d'], location)
            shutil.copy(motif['2d'], location)
            self.logger.debug("Copying motif 2d %s %s to %s",
                              motif['motif_id'], motif['2d'], web_location)
            shutil.copy(motif['2d'], web_location)
            self.logger.debug("To be copying motif mat file %s %s to %s",
                              motif['motif_id'], mat_file, mat_out)
            shutil.copy(mat_file, mat_out) 

        return None
