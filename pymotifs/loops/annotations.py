import csv

from pymotifs import core
from pymotifs import models as mod
from sqlalchemy import desc

class Loader(core.SimpleLoader):

    merge_data = True
    allow_no_data = True
    mark = False #Add this line so that each time this stage is run, it will re-process items in to_process()

    """
    @property
    def table(self):
        return mod.LoopAnnotations
    """

    def to_process(self, pdbs, **kwrags):
        """ return the different loop types to process, see if that works
        """
        return ["IL"]

    """
    def remove(self, *args, **kwargs):
        Does not actually remove from the DB. We always want the loop ids to
        be consistent so we do not automatically remove loops.
        """

    def query(self, session, pdb):
        return session.query(mod.LoopAnnotations)

    def _get_loop_annotations(self):
        """ Gets access to the spreadsheet, reads values, and saves data in DB friendly format.
        """
        #Create dictionary for loop_annotations and motif_annotations with loop_id and motif_id as their respective keys and values consisting of a 3-element set: annotation_1, annotation_2, and author.

        loop_annotations = {}
        motif_annotations = {}

        #Read data from csv file and save each row as a list of strings

        print("annotations.py reading loop annotations")

        with open("/usr/local/pipeline/hub-core/IL_annotations.csv", "r") as annotations:
            entries = csv.reader(annotations, delimiter = ',')
            for row in entries:
                if row[1] != "":
                    loop_annotations[row[1]] = [row[2], row[3], row[4]]
                else:
                    motif_annotations[row[0]] = [row[2], row[3], row[4]]

        #Save each loop annnotation to loop_data in DB-friendly format
        if not loop_annotations:
            self.logger.error("No loop annotations found.")
        else:
            loop_data = []
            for loop in loop_annotations.keys():
                loop_data.append(mod.LoopAnnotations(
                    loop_id = loop,
                    annotation_1 = loop_annotations[loop][0],
                    annotation_2 = loop_annotations[loop][1],
                    author = loop_annotations[loop][2]))

        if not motif_annotations:
            self.logger.error("No motif annotations found.")
        else:
            with self.session() as session:
                current_ml_release = session.query(mod.MlReleases.ml_release_id).\
                    order_by(desc(mod.MlReleases.date)).\
                    limit(1)

            print(current_ml_release)

            #For each motif_annotation, get the associated loop_ids. Exclude loop_ids with their own annotations. Create motif_annotation in DB-friendly format for the remaining loop_ids.
            for motif in motif_annotations.keys():
                with self.session() as session:
                    query = session.query(mod.MlLoops.loop_id).\
                         filter(mod.MlLoops.motif_id == motif).\
                         filter(mod.MlLoops.ml_release_id == current_ml_release)
                    loop_ids = [r.loop_id for r in query]

                loop_ids = set(loop_ids).difference(loop_annotations.keys())
                for loop in loop_ids:
                    loop_data.append(mod.LoopAnnotations(
                        loop_id = loop,
                        annotation_1 = motif_annotations[motif][0],
                        annotation_2 = motif_annotations[motif][1],
                        author = motif_annotations[motif][2]))
        return loop_data

    def data(self, pdb, **kwargs):
        """ Get the loop annotations from the spreadsheet.
        Return a list of all loop annotations. """

        print("annotations.py running data method")

        #Get loop annotations
        data = []
        data = self._get_loop_annotations()
#Return data to be written to the database.
        if not data:
            self.logger.error("No loop annotations found.")
        else:
            return data
