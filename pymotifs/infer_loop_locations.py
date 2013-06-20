import sys
import logging

from sqlalchemy import distinct

from models import session
from models import Feature
from models import AllLoops
from models import FeatureType
from models import LoopLocation
from models import FeatureNucleotide
from models import LoopLocationAnnotation
from MotifAtlasBaseClass import MotifAtlasBaseClass


class InferLocations(MotifAtlasBaseClass):
    """A class to infer the location annotations for loops. The basic idea is
    that we get all loops from some structure and all annotated features and
    then seeing what nucleotides of the loop overlap with annotated features
    and then using the coresponding feature types to build up loop location
    annotations.
    """

    def load(self, pdbs):
        """This will load all the loop location annotations for all given pdbs.
        If no pdbs are given then this will find all pdbs with feature
        annotations and then mark those. This will always recompute this for
        all pdbs ignoring the data in PdbAnalysisStatus. This is because we
        don't have a good way of noticing when we have changed the location
        names and the like. This should be changed in the future so we do not
        always recompute.
        """

        pdbs = pdbs or self.__find_annotated_pdbs__()
        logging.info("Infering locations for %s pdbs", len(pdbs))

        for pdb in pdbs:
            logging.info("Infering location for %s", pdb)
            self.process_pdb(pdb)

    def process_pdb(self, pdb):
        """Infer the location annotations for a single pdb.
        """

        loops = self.loops(pdb)
        ids = [loop.id for loop in loops]
        query = session.query(LoopLocationAnnotation).\
            filter(LoopLocationAnnotation.loop_id.in_(ids))

        if query.count() != 0:
            logging.info("Removing previous annotations for %s", pdb)
            query.delete(synchronize_session=False)

        for loop in loops:
            logging.info("Processing loop: %s", loop.id)

            features = self.loop_features(loop)
            locations = self.loop_locations(loop, features)

            annotations = self.generate_annotations(loop, locations)
            if len(annotations) == 0:
                logging.error("Did not find any annotaitons for %s", loop.id)
            else:
                logging.info("Found %s annotations for %s", len(annotations),
                             loop.id)

            logging.info("Storing all annotations")

            try:
                for annotation in annotations:
                    session.add(annotation)
                session.commit()
            except:
                logging.warn("Failed to store all annotations for %s",
                             loop.name)
                return False

            logging.info("Stored locations for %s", pdb)
        self.mark_pdb_as_analyzed(pdb, 'loop_annotation')

    def loops(self, pdb):
        """Get all loops in a pdb.
        """
        return session.query(AllLoops).\
            filter(AllLoops.pdb == pdb)

    def loop_locations(self, loop, types):
        """For each given feature type we generate a loop location and then
        store it if the generated name is new.
        """

        locations = []
        for feature_type in types:
            location = self.loop_location(loop, feature_type)
            if not location:
                continue

            query = session.query(LoopLocation).\
                filter(LoopLocation.name == location.name)

            count = query.count()
            if count == 0:
                logging.info("Add new loop location: %s", location.name)
                session.add(location)
                locations.append(location)
            elif count == 1:
                locations.append(query.first())
            else:
                logging.error("Mulitple loops with same name: %s.",
                              location.name)

        session.commit()
        return locations

    def loop_location(self, loop, feature_type):
        """Create the LoopLocation for the given loop and feature type. In the
        case where the loop is of an unknown type, ie not IL, HL or J3 we log
        an error and return None.
        """

        location = feature_type.loop_location(loop)

        if not location:
            logging.error("Unknown loop type %s for %s", loop.type, loop.id)
            return None

        return location

    def loop_features(self, loop):
        """Get all feature type for the given loop. This is done by finding all
        features that the loops nucleotides match and then finding all feature
        types of the matched features.
        """

        query = session.query(FeatureType).\
            join(Feature, Feature.type_id == FeatureType.id).\
            join(FeatureNucleotide,
                 FeatureNucleotide.feature_id == Feature.id).\
            filter(FeatureNucleotide.unit_id.in_(loop.nucleotides()))

        return list(set(query))

    def generate_annotations(self, loop, locations):
        """Gernate the LoopLocationAnnotations for each location using the
        given loop.
        """

        annotations = []
        for location in locations:
            annotation = LoopLocationAnnotation(loop_id=loop.id,
                                                loop_location_id=location.id)
            annotations.append(annotation)
        return annotations

    def __find_annotated_pdbs__(self):
        logging.info("Finding all pdbs with annotated locations")
        query = session.query(distinct(Feature.pdb))
        return [feature[0] for feature in query]


def main(pdbs):
    loader = InferLocations()
    loader.start_logging()
    loader.load(pdbs)

if __name__ == '__main__':
    main(sys.argv[1:])
